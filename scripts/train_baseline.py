"""
Train models to predict Toxoplasma gondii IgG seropositivity (NHANES 2013-2014 adults).

Input: data/processed/nhanes_adult_toxo_clean.csv (built by scripts/cleaning2.py)

Steps:
  1. Load data, define feature roles, stratified 80/20 train/test split
  2. Preprocessing pipeline (imputation, missing flags, one-hot, scaling)
  3. Model comparison: baselines + 3 models x 3 feature sets, repeated CV on
     the training set only, with hypothesis tests H1-H4
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.compose import ColumnTransformer
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold, cross_validate, train_test_split
from sklearn.pipeline import Pipeline, make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from xgboost import XGBClassifier

DATA_CSV = os.path.join('data', 'processed', 'nhanes_adult_toxo_clean.csv')
TARGET = 'target_toxo_pos'
RANDOM_STATE = 42
TEST_SIZE = 0.20

# ---- 1. Load, feature roles, split ------------------------------------------
# Never used as features: ID, survey design, and DMDBORN4 (duplicated by
# yrs_in_us_grp, whose level 0 = US-born)
NOT_FEATURES = ['SEQN', 'WTMEC2YR', 'SDMVPSU', 'SDMVSTRA', 'DMDBORN4']

#Nominal codes
CATEGORICAL = ['RIAGENDR', 'RIDRETH3', 'yrs_in_us_grp', 'SLQ050', 'smoking_status',
               'toxocara_pos']
#Continuous or ordered codes
NUMERIC = ['RIDAGEYR', 'DMDEDUC2', 'INDFMPIR', 'HSD010', 'SLD010H',
           'drinks_per_day', 'phq9_total_score']

DEMOGRAPHICS = ['RIDAGEYR', 'RIAGENDR', 'RIDRETH3', 'yrs_in_us_grp', 'DMDEDUC2', 'INDFMPIR']
BEHAVIOUR = ['HSD010', 'SLD010H', 'SLQ050', 'smoking_status', 'drinks_per_day', 'phq9_total_score']

# Nested sets: each one adds a single group to the previous one
FEATURE_SETS = {
    'age_only': ['RIDAGEYR'],                                       # simplest real baseline
    'A_demographics': DEMOGRAPHICS,
    'B_demo_toxocara': DEMOGRAPHICS + ['toxocara_pos'],
    'C_all': DEMOGRAPHICS + ['toxocara_pos'] + BEHAVIOUR,
}
assert set(FEATURE_SETS['C_all']) == set(CATEGORICAL + NUMERIC), "C_all must use every feature"

df = pd.read_csv(DATA_CSV)
unassigned = set(df.columns) - set(NOT_FEATURES) - set(CATEGORICAL) - set(NUMERIC) - {TARGET}
assert not unassigned, f"Columns without a role: {unassigned}"

X = df[CATEGORICAL + NUMERIC]
y = df[TARGET]
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=TEST_SIZE, stratify=y, random_state=RANDOM_STATE)

print(f"Loaded {len(df)} rows, {X.shape[1]} candidate features")
print(f"Train: {len(X_train)} rows, {y_train.mean():.1%} positive")
print(f"Test:  {len(X_test)} rows, {y_test.mean():.1%} positive  (held out until the final step)")
print("\nMissing values in the training set (imputed in the pipeline):")
print(X_train.isna().sum()[lambda s: s > 0].to_string())

# ---- 2. Preprocessing ---------------------------------------------------------
# Numeric columns with substantial missingness also get a 0/1 "was missing"
# column (missingness is not random here). DMDEDUC2/SLD010H have only a few
# missing values, so a flag would be a near-constant column -> median only.
NUMERIC_FLAG_MISSING = ['INDFMPIR', 'HSD010', 'drinks_per_day', 'phq9_total_score']


def make_preprocessor(features):
    """Fit-inside-the-pipeline preprocessing for the given feature list."""
    cat = [c for c in CATEGORICAL if c in features]
    num_flag = [c for c in NUMERIC if c in features and c in NUMERIC_FLAG_MISSING]
    num_plain = [c for c in NUMERIC if c in features and c not in NUMERIC_FLAG_MISSING]
    steps = []
    if cat:
        steps.append(('cat', Pipeline([
            ('impute', SimpleImputer(strategy='most_frequent')),
            ('onehot', OneHotEncoder(drop='if_binary', handle_unknown='ignore',
                                     sparse_output=False)),
        ]), cat))
    if num_plain:
        steps.append(('num', Pipeline([
            ('impute', SimpleImputer(strategy='median')),
            ('scale', StandardScaler()),
        ]), num_plain))
    if num_flag:
        steps.append(('numflag', Pipeline([
            ('impute', SimpleImputer(strategy='median', add_indicator=True)),
            ('scale', StandardScaler()),
        ]), num_flag))
    return ColumnTransformer(steps)


# Verify the preprocessor on the training set only (the test set stays untouched)
print("\n--- Step 2 check: preprocessing fitted on the training set ---")
for name, cols in FEATURE_SETS.items():
    pre = make_preprocessor(cols).set_output(transform='pandas')
    Xt = pre.fit_transform(X_train[cols])
    assert Xt.isna().sum().sum() == 0, "NaNs left after preprocessing"
    assert len(Xt) == len(X_train)
    print(f"\n[{name}] {len(cols)} input columns -> {Xt.shape[1]} model columns:")
    print("  " + ", ".join(Xt.columns))

    # missing flags: column sum must equal the raw missing count
    flag_cols = [c for c in Xt.columns if 'missingindicator' in c]
    for fc in flag_cols:
        raw = fc.split('missingindicator_')[1]
        n_flag = int((Xt[fc] > Xt[fc].min()).sum())
        assert n_flag == X_train[raw].isna().sum(), f"flag mismatch for {raw}"
    print(f"  missing flags match raw missing counts: {len(flag_cols)} checked")

    # scaled columns: mean ~0, std ~1
    scaled = [c for c in Xt.columns if not c.startswith('cat__')]
    print(f"  scaled columns mean range [{Xt[scaled].mean().min():.3f}, {Xt[scaled].mean().max():.3f}], "
          f"std range [{Xt[scaled].std(ddof=0).min():.3f}, {Xt[scaled].std(ddof=0).max():.3f}]")

    # one-hot columns: each row must have exactly one level per multi-level variable
    for var in [c for c in CATEGORICAL if c in cols and X_train[c].nunique() > 2]:
        block = Xt[[c for c in Xt.columns if c.startswith(f'cat__{var}_')]]
        assert (block.sum(axis=1) == 1).all(), f"one-hot rows not exactly one level for {var}"
    print("  one-hot blocks: exactly one level per row")

# Medians used for imputation (learned from training rows only)
pre_all = make_preprocessor(FEATURE_SETS['C_all']).fit(X_train)
medians = dict(zip(pre_all.named_transformers_['num'].feature_names_in_,
                   pre_all.named_transformers_['num'].named_steps['impute'].statistics_))
medians.update(zip(pre_all.named_transformers_['numflag'].feature_names_in_,
                   pre_all.named_transformers_['numflag'].named_steps['impute'].statistics_))
print("\nImputation medians (train only):", {k: round(float(v), 2) for k, v in medians.items()})

# ---- 3. Model comparison (cross-validation on the training set only) ---------
# One fixed list of splits -> every model is scored on exactly the same 15 folds
CV_SPLITS = list(RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=RANDOM_STATE)
                 .split(X_train, y_train))

MODELS = {
    'Dummy': DummyClassifier(strategy='prior'),                       # same score for everyone
    'LogReg': LogisticRegression(C=1.0, max_iter=5000),
    'RandomForest': RandomForestClassifier(n_estimators=500, min_samples_leaf=10,
                                           n_jobs=-1, random_state=RANDOM_STATE),
    'GradBoost': HistGradientBoostingClassifier(learning_rate=0.05, max_depth=3, max_iter=300,
                                                random_state=RANDOM_STATE),
    # Conservative settings mirroring GradBoost (shallow trees, small learning rate)
    'XGBoost': XGBClassifier(n_estimators=300, learning_rate=0.05, max_depth=3,
                             subsample=0.8, colsample_bytree=0.8, eval_metric='logloss',
                             n_jobs=-1, random_state=RANDOM_STATE),
}
REAL_MODELS = ['LogReg', 'RandomForest', 'GradBoost', 'XGBoost']
EXPERIMENTS = [('Dummy', 'age_only'), ('LogReg', 'age_only')] + [
    (m, s) for s in ['A_demographics', 'B_demo_toxocara', 'C_all'] for m in REAL_MODELS]

print("\n--- Step 3 check: cross-validation folds ---")
assert set(X_train.index).isdisjoint(X_test.index), "test rows leaked into training data"
val_counts = np.zeros(len(X_train), dtype=int)
for _, val_idx in CV_SPLITS:
    val_counts[val_idx] += 1
assert (val_counts == 3).all(), "every training row must be validated exactly once per repeat"
fold_pos = [y_train.iloc[v].mean() for _, v in CV_SPLITS]
print(f"  {len(CV_SPLITS)} folds; every training row validated exactly 3 times; "
      f"positive share per validation fold {min(fold_pos):.1%}-{max(fold_pos):.1%}; test set not used")

results = {}
for model_name, set_name in EXPERIMENTS:
    cols = FEATURE_SETS[set_name]
    pipe = make_pipeline(make_preprocessor(cols), MODELS[model_name])
    cv = cross_validate(pipe, X_train[cols], y_train, cv=CV_SPLITS,
                        scoring={'roc_auc': 'roc_auc', 'pr_auc': 'average_precision',
                                 'brier': 'neg_brier_score'})
    results[(model_name, set_name)] = {'roc_auc': cv['test_roc_auc'], 'pr_auc': cv['test_pr_auc'],
                                       'brier': -cv['test_brier']}

dummy = results[('Dummy', 'age_only')]
assert np.allclose(dummy['roc_auc'], 0.5), "dummy must score exactly 0.5 ROC-AUC"
expected_brier = y_train.mean() * (1 - y_train.mean())
print(f"  dummy ROC-AUC = 0.500 in all folds; dummy Brier {dummy['brier'].mean():.4f} "
      f"(theory p(1-p) = {expected_brier:.4f})")

print("\nCV results (mean +/- sd over 15 folds; Brier: lower is better)")
table = pd.DataFrame([{
    'model': m, 'features': s,
    'ROC-AUC': f"{r['roc_auc'].mean():.3f} +/- {r['roc_auc'].std():.3f}",
    'PR-AUC': f"{r['pr_auc'].mean():.3f} +/- {r['pr_auc'].std():.3f}",
    'Brier': f"{r['brier'].mean():.4f}",
} for (m, s), r in results.items()])
print(table.to_string(index=False))


def compare(base, other, label):
    """Paired ROC-AUC difference (other - base) on the same folds, with the
    corrected repeated-CV t-test (Nadeau & Bengio), since CV folds overlap."""
    d = results[other]['roc_auc'] - results[base]['roc_auc']
    k, test_train_ratio = len(d), 0.2 / 0.8
    t = d.mean() / np.sqrt((1 / k + test_train_ratio) * d.var(ddof=1))
    p = 2 * stats.t.sf(abs(t), df=k - 1)
    meaningful = d.mean() >= 0.01 and p < 0.05
    print(f"  {label:44} diff {d.mean():+.4f}  better in {(d > 0).sum():2}/{k} folds  "
          f"p={p:.3f}  -> {'MEANINGFUL' if meaningful else 'not meaningful'}")


print("\nHypothesis tests (meaningful = ROC-AUC gain >= 0.01 AND p < 0.05)")
print("H1: real models beat the baselines (feature set C)")
for m in REAL_MODELS:
    compare(('Dummy', 'age_only'), (m, 'C_all'), f"{m} vs Dummy")
    compare(('LogReg', 'age_only'), (m, 'C_all'), f"{m} vs age-only LogReg")
print("H2: complex models beat logistic regression (expected: no)")
for s in ['A_demographics', 'B_demo_toxocara', 'C_all']:
    for m in REAL_MODELS[1:]:
        compare(('LogReg', s), (m, s), f"{m} vs LogReg [{s}]")
print("H3: Toxocara adds to demographics (B vs A; expected: small gain)")
for m in REAL_MODELS:
    compare((m, 'A_demographics'), (m, 'B_demo_toxocara'), f"{m}: B vs A")
print("H4: behaviour/health adds to B (C vs B; expected: no)")
for m in REAL_MODELS:
    compare((m, 'B_demo_toxocara'), (m, 'C_all'), f"{m}: C vs B")
