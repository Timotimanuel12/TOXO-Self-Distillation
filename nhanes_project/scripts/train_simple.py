"""
Simple training script: predict Toxoplasma infection (NHANES 2013-2014 adults).

Input: data/processed/nhanes_adult_toxo_clean_v2.csv (made by scripts/cleaning3.py)
(v2 = v1 + the lifestyle features from the IEEE HIV-cohort paper; same people, same order)

Steps:
  1. Load the clean data and split it: 80% training, 20% test
  2. Prepare the data: fill gaps, add "was missing" flags, one-hot encode, scale
  3. Compare models with cross-validation (training set only)
  4. Test the hypotheses H1-H5
"""
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.compose import ColumnTransformer
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold, cross_validate, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from xgboost import XGBClassifier


# =====================================================================
# STEP 1: Load the data and split it
# =====================================================================
df = pd.read_csv('data/processed/nhanes_adult_toxo_clean_v2.csv')
# Lifetime partners runs from 0 to 1000, so use log(1 + x) to shrink the long tail
df['log_lifetime_partners'] = np.log1p(df['lifetime_partners'])

# Category columns: the numbers are labels (e.g. race 2 = Other Hispanic)
CATEGORICAL = ['RIAGENDR', 'RIDRETH3', 'yrs_in_us_grp', 'SLQ050', 'smoking_status', 'toxocara_pos']

# Lifestyle features from the IEEE HIV-cohort paper (sexual behaviour, drug use, oral health, HIV).
# inject_count_grp is left out: it is almost an exact copy of inject_drugs_ever.
PAPER = ['anal_sex_ever', 'same_sex_partner_ever', 'log_lifetime_partners', 'sex_freq_year',
         'condomless_sex_freq', 'oral_sex_no_protection', 'age_first_sex', 'std_ever',
         'inject_drugs_ever', 'hard_drug_ever',
         'missing_teeth', 'loose_teeth', 'mouth_aching_freq', 'hiv_pos']

# Number columns: the numbers are real amounts (e.g. age in years)
NUMERIC = ['RIDAGEYR', 'DMDEDUC2', 'INDFMPIR', 'HSD010', 'SLD010H', 'drinks_per_day', 'phq9_total_score'] + PAPER

# Number columns with many empty cells: these also get a "was missing" flag
# (the PAPER questions were only asked in some age groups, so they all get one)
FLAG_MISSING = ['INDFMPIR', 'HSD010', 'drinks_per_day', 'phq9_total_score'] + PAPER

# Feature sets: each one adds one group to the previous one
DEMOGRAPHICS = ['RIDAGEYR', 'RIAGENDR', 'RIDRETH3', 'yrs_in_us_grp', 'DMDEDUC2', 'INDFMPIR']
BEHAVIOUR = ['HSD010', 'SLD010H', 'SLQ050', 'smoking_status', 'drinks_per_day', 'phq9_total_score']
FEATURE_SETS = {
    'age_only': ['RIDAGEYR'],
    'A_demographics': DEMOGRAPHICS,
    'B_demo_toxocara': DEMOGRAPHICS + ['toxocara_pos'],
    'C_all': DEMOGRAPHICS + ['toxocara_pos'] + BEHAVIOUR,
    'D_paper': DEMOGRAPHICS + ['toxocara_pos'] + BEHAVIOUR + PAPER,
}

# Features (X) and the answer we want to predict (y)
# (SEQN, the survey columns and DMDBORN4 are never used as features)
X = df[CATEGORICAL + NUMERIC]
y = df['target_toxo_pos']

# 80% training, 20% test. stratify=y keeps 15% positive in both parts.
# random_state=42 makes the split the same every time.
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.20, stratify=y, random_state=42)

print('Training set:', len(X_train), 'people,', y_train.sum(), 'positive')
print('Test set:    ', len(X_test), 'people,', y_test.sum(), 'positive (not used until the final step)')


# =====================================================================
# STEP 2: Data preparation recipe
# =====================================================================
def make_preprocessor(features):
    """Build the preparation recipe for a list of feature columns."""

    # Sort the chosen features into three groups
    cat = []          # category columns
    num_flag = []     # number columns with many gaps (get a "was missing" flag)
    num_plain = []    # number columns with few or no gaps
    for col in CATEGORICAL:
        if col in features:
            cat.append(col)
    for col in NUMERIC:
        if col in features:
            if col in FLAG_MISSING:
                num_flag.append(col)
            else:
                num_plain.append(col)

    # Category columns: fill gaps with the most common value, then one-hot encode
    cat_recipe = Pipeline([
        ('fill', SimpleImputer(strategy='most_frequent')),
        ('onehot', OneHotEncoder(drop='if_binary', handle_unknown='ignore', sparse_output=False)),
    ])
    # Number columns with few gaps: fill gaps with the median, then scale
    plain_recipe = Pipeline([
        ('fill', SimpleImputer(strategy='median')),
        ('scale', StandardScaler()),
    ])
    # Number columns with many gaps: fill with the median + add a "was missing" flag, then scale
    flag_recipe = Pipeline([
        ('fill', SimpleImputer(strategy='median', add_indicator=True)),
        ('scale', StandardScaler()),
    ])

    # Send each group of columns to its recipe (skip groups with no columns)
    groups = []
    if len(cat) > 0:
        groups.append(('cat', cat_recipe, cat))
    if len(num_plain) > 0:
        groups.append(('num', plain_recipe, num_plain))
    if len(num_flag) > 0:
        groups.append(('numflag', flag_recipe, num_flag))
    return ColumnTransformer(groups)


# =====================================================================
# STEP 3: Compare models with cross-validation (training set only)
# =====================================================================
MODELS = {
    'Dummy': DummyClassifier(strategy='prior'),   # gives everyone the same 15% risk
    'LogReg': LogisticRegression(C=1.0, max_iter=5000),
    'RandomForest': RandomForestClassifier(n_estimators=500, min_samples_leaf=10,
                                           n_jobs=-1, random_state=42),
    'GradBoost': HistGradientBoostingClassifier(learning_rate=0.05, max_depth=3, max_iter=300,
                                                random_state=42),
    'XGBoost': XGBClassifier(n_estimators=300, learning_rate=0.05, max_depth=3,
                             subsample=0.8, colsample_bytree=0.8, eval_metric='logloss',
                             n_jobs=-1, random_state=42),
}
REAL_MODELS = ['LogReg', 'RandomForest', 'GradBoost', 'XGBoost']

# Which model is trained on which feature set
experiments = [('Dummy', 'age_only'), ('LogReg', 'age_only')]
for set_name in ['A_demographics', 'B_demo_toxocara', 'C_all', 'D_paper']:
    for model_name in REAL_MODELS:
        experiments.append((model_name, set_name))

# 5 groups x 3 repeats = 15 practice tests. Made once, so every model gets the same tests.
cv_splits = list(RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=42)
                 .split(X_train, y_train))

print('\nCV results (mean +/- sd over 15 practice tests; Brier: lower is better)')
results = {}
for model_name, set_name in experiments:
    cols = FEATURE_SETS[set_name]
    pipe = Pipeline([
        ('prepare', make_preprocessor(cols)),
        ('model', MODELS[model_name]),
    ])
    # Train on 4/5 and test on 1/5, 15 times
    scores = cross_validate(pipe, X_train[cols], y_train, cv=cv_splits,
                            scoring=['roc_auc', 'average_precision', 'neg_brier_score'])
    roc = scores['test_roc_auc']
    pr = scores['test_average_precision']
    brier = -scores['test_neg_brier_score']
    results[(model_name, set_name)] = roc

    print(f'{model_name:13} {set_name:16} ROC-AUC {roc.mean():.3f} +/- {roc.std():.3f}   '
          f'PR-AUC {pr.mean():.3f} +/- {pr.std():.3f}   Brier {brier.mean():.4f}')


# =====================================================================
# STEP 4: Hypothesis tests
# =====================================================================
def compare(base, other, label):
    """Is 'other' better than 'base'? Compares ROC-AUC on the same 15 practice tests.

    Meaningful = at least +0.01 ROC-AUC and p < 0.05. The p-value uses the corrected
    t-test (Nadeau & Bengio), because the practice tests share training people.
    """
    diff = results[other] - results[base]
    k = len(diff)
    t = diff.mean() / np.sqrt((1 / k + 0.2 / 0.8) * diff.var(ddof=1))
    p = 2 * stats.t.sf(abs(t), df=k - 1)
    better_in = (diff > 0).sum()

    if diff.mean() >= 0.01 and p < 0.05:
        verdict = 'MEANINGFUL'
    else:
        verdict = 'not meaningful'
    print(f'  {label:44} diff {diff.mean():+.4f}  better in {better_in:2}/{k} tests  '
          f'p={p:.3f}  -> {verdict}')


print('\nHypothesis tests (meaningful = ROC-AUC gain >= 0.01 AND p < 0.05)')

print('H1: real models beat the baselines (feature set C)')
for m in REAL_MODELS:
    compare(('Dummy', 'age_only'), (m, 'C_all'), m + ' vs Dummy')
    compare(('LogReg', 'age_only'), (m, 'C_all'), m + ' vs age-only LogReg')

print('H2: complex models beat logistic regression (expected: no)')
for s in ['A_demographics', 'B_demo_toxocara', 'C_all']:
    for m in ['RandomForest', 'GradBoost', 'XGBoost']:
        compare(('LogReg', s), (m, s), m + ' vs LogReg [' + s + ']')

print('H3: Toxocara adds to demographics (B vs A; expected: small gain)')
for m in REAL_MODELS:
    compare((m, 'A_demographics'), (m, 'B_demo_toxocara'), m + ': B vs A')

print('H4: behaviour/health adds to B (C vs B; expected: no)')
for m in REAL_MODELS:
    compare((m, 'B_demo_toxocara'), (m, 'C_all'), m + ': C vs B')

print('H5: lifestyle features from the HIV-cohort paper add to C (D vs C; expected: no)')
for m in REAL_MODELS:
    compare((m, 'C_all'), (m, 'D_paper'), m + ': D vs C')
