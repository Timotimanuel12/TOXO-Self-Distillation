"""
Final evaluation: the one and only use of the 20% test set (RQ4 and RQ5).

Input: data/processed/nhanes_adult_toxo_clean_v2.csv (same 80/20 split as scripts/train_simple.py)

Final model: logistic regression on Set A (demographics). In train_simple.py no feature group
(H3 Toxocara, H4 behaviour/health, H5 HIV-paper features) and no complex model (H2) gave a
meaningful gain, so the simpler option is kept. It needs questionnaire answers only (no blood test).

Steps:
  1. Same 80/20 split as train_simple.py
  2. Tune C (training set only, same 15 practice tests)
  3. Choose a risk cut-off from out-of-fold training predictions (Youden's J)
  4. Fit on the whole training set and score the test set ONCE
  5. Discrimination: ROC-AUC, PR-AUC, Brier with bootstrap 95% CIs, against both baselines
  6. Calibration: intercept, slope, risk-group table and plot
  7. Cut-off metrics and decision curve
  8. RQ4: adjusted odds ratios (training set, all Set C variables)
  9. RQ4: SHAP values for the final model on the test set
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
from scipy import stats
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score, roc_curve
from sklearn.model_selection import (RepeatedStratifiedKFold, StratifiedKFold, cross_val_predict,
                                     cross_validate, train_test_split)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

FIG_DIR = os.path.join('outputs', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

# Chart style: light surface, recessive grid and axes, text in ink colours (not series colours)
INK, INK2, GRID, SURFACE = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
BLUE, ORANGE = '#2a78d6', '#eb6834'
plt.rcParams.update({
    'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE,
    'axes.edgecolor': INK2, 'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2,
    'axes.titlecolor': INK, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.8, 'font.size': 10,
})


# =====================================================================
# STEP 1: Load the data and make the same split as train_simple.py
# =====================================================================
df = pd.read_csv('data/processed/nhanes_adult_toxo_clean_v2.csv')
y = df['target_toxo_pos']

# Same arguments as train_simple.py, so the same 4,004 / 1,002 people
train, test = train_test_split(df, test_size=0.20, stratify=y, random_state=42)
y_train, y_test = train['target_toxo_pos'], test['target_toxo_pos']
print('Training set:', len(train), 'people,', y_train.sum(), 'positive')
print('Test set:    ', len(test), 'people,', y_test.sum(), 'positive')

# Set A (demographics) = the final model's features
CATEGORICAL = ['RIAGENDR', 'RIDRETH3', 'yrs_in_us_grp']
NUMERIC = ['RIDAGEYR', 'DMDEDUC2']
FLAG_MISSING = ['INDFMPIR']
FINAL = CATEGORICAL + NUMERIC + FLAG_MISSING

LABELS = {'RIDAGEYR': 'Age', 'RIDRETH3': 'Race/ethnicity', 'yrs_in_us_grp': 'Birthplace / years in US',
          'DMDEDUC2': 'Education', 'INDFMPIR': 'Income-to-poverty ratio', 'RIAGENDR': 'Sex'}


def make_pipeline(C=1.0):
    """Same preparation recipe as train_simple.py, for the Set A columns."""
    prepare = ColumnTransformer([
        ('cat', Pipeline([('fill', SimpleImputer(strategy='most_frequent')),
                          ('onehot', OneHotEncoder(drop='if_binary', handle_unknown='ignore',
                                                   sparse_output=False))]), CATEGORICAL),
        ('num', Pipeline([('fill', SimpleImputer(strategy='median')),
                          ('scale', StandardScaler())]), NUMERIC),
        ('numflag', Pipeline([('fill', SimpleImputer(strategy='median', add_indicator=True)),
                              ('scale', StandardScaler())]), FLAG_MISSING),
    ])
    return Pipeline([('prepare', prepare), ('model', LogisticRegression(C=C, max_iter=5000))])


def corrected_p(diff):
    """Nadeau-Bengio corrected t-test, as in train_simple.py."""
    k = len(diff)
    t = diff.mean() / np.sqrt((1 / k + 0.2 / 0.8) * diff.var(ddof=1))
    return 2 * stats.t.sf(abs(t), df=k - 1)


# =====================================================================
# STEP 2: Tune C (how strongly the model is held back), training set only
# =====================================================================
# Keep C = 1 unless another value is meaningfully better (same rule as everywhere else)
cv_splits = list(RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=42)
                 .split(train[FINAL], y_train))
print('\nTuning C (CV ROC-AUC, training set only)')
tuning = {}
for C in [0.001, 0.01, 0.1, 1.0, 10.0, 100.0]:
    tuning[C] = cross_validate(make_pipeline(C), train[FINAL], y_train, cv=cv_splits,
                               scoring='roc_auc')['test_score']
for C, scores in tuning.items():
    diff = scores - tuning[1.0]
    p = corrected_p(diff) if C != 1.0 else float('nan')
    print(f'  C={C:<7} ROC-AUC {scores.mean():.4f}   vs C=1: {diff.mean():+.4f}  p={p:.3f}')

best_C = 1.0
for C, scores in tuning.items():
    diff = scores - tuning[1.0]
    if C != 1.0 and diff.mean() >= 0.01 and corrected_p(diff) < 0.05:
        best_C = C
print(f'  -> chosen C = {best_C}')


# =====================================================================
# STEP 3: Choose a risk cut-off from out-of-fold training predictions
# =====================================================================
# Each training person gets a prediction from a model that did not see them.
# The cut-off maximises Youden's J (sensitivity + specificity - 1). The test set is not used.
oof = cross_val_predict(make_pipeline(best_C), train[FINAL], y_train,
                        cv=StratifiedKFold(5, shuffle=True, random_state=42),
                        method='predict_proba')[:, 1]
fpr, tpr, thresholds = roc_curve(y_train, oof)
cutoff = thresholds[np.argmax(tpr - fpr)]
print(f'\nRisk cut-off chosen on training data (Youden): {cutoff:.3f}')


# =====================================================================
# STEP 4: Fit on the whole training set, then score the test set ONCE
# =====================================================================
final = make_pipeline(best_C).fit(train[FINAL], y_train)
p_test = final.predict_proba(test[FINAL])[:, 1]
yt = y_test.to_numpy()

# Baselines, also fitted on the training set only
p_dummy = np.full(len(test), y_train.mean())
age_only = Pipeline([('scale', StandardScaler()), ('model', LogisticRegression(max_iter=5000))])
p_age = age_only.fit(train[['RIDAGEYR']], y_train).predict_proba(test[['RIDAGEYR']])[:, 1]


# =====================================================================
# STEP 5: Discrimination with bootstrap 95% confidence intervals
# =====================================================================
rng = np.random.default_rng(42)
boot_idx = [rng.integers(0, len(yt), len(yt)) for _ in range(2000)]
boot_idx = [i for i in boot_idx if 0 < yt[i].sum() < len(i)]   # each resample needs both classes


def with_ci(metric, p):
    """Metric on the test set + 2.5th and 97.5th percentiles over 2,000 bootstrap resamples."""
    vals = [metric(yt[i], p[i]) for i in boot_idx]
    return metric(yt, p), np.percentile(vals, 2.5), np.percentile(vals, 97.5)


print('\nTest set results (95% CI from 2,000 bootstrap resamples)')
for name, p in [('Dummy (prevalence)', p_dummy), ('Age-only LR', p_age), ('Final model (LR, Set A)', p_test)]:
    roc = with_ci(roc_auc_score, p) if name != 'Dummy (prevalence)' else (0.5, 0.5, 0.5)
    pr = with_ci(average_precision_score, p)
    br = with_ci(brier_score_loss, p)
    print(f'  {name:24} ROC-AUC {roc[0]:.3f} ({roc[1]:.3f}-{roc[2]:.3f})   '
          f'PR-AUC {pr[0]:.3f} ({pr[1]:.3f}-{pr[2]:.3f})   Brier {br[0]:.4f} ({br[1]:.4f}-{br[2]:.4f})')

# Is the final model better than age alone on the test set? (paired bootstrap)
diffs = [roc_auc_score(yt[i], p_test[i]) - roc_auc_score(yt[i], p_age[i]) for i in boot_idx]
print(f'  ROC-AUC gain over age-only: {roc_auc_score(yt, p_test) - roc_auc_score(yt, p_age):+.3f} '
      f'(95% CI {np.percentile(diffs, 2.5):+.3f} to {np.percentile(diffs, 97.5):+.3f})')
print(f'  Cross-validated ROC-AUC on the training set was {tuning[best_C].mean():.3f}')


# =====================================================================
# STEP 6: Calibration (do predicted risks match observed rates?)
# =====================================================================
logit_p = np.log(p_test / (1 - p_test)).reshape(-1, 1)
recal = LogisticRegression(penalty=None).fit(logit_p, yt)
print('\nCalibration on the test set')
print(f'  Mean predicted risk {p_test.mean():.3f}   observed rate {yt.mean():.3f}')
print(f'  Calibration slope {recal.coef_[0][0]:.2f} (ideal 1)   intercept {recal.intercept_[0]:+.2f} (ideal 0)')

groups = pd.qcut(p_test, 10, labels=False, duplicates='drop')
calib = (pd.DataFrame({'group': groups, 'p': p_test, 'y': yt})
         .groupby('group').agg(n=('y', 'size'), predicted=('p', 'mean'), observed=('y', 'mean')))
print('  Risk groups (tenths of predicted risk):')
print(calib.round(3).to_string())

fig, ax = plt.subplots(figsize=(5.2, 4.6))
top = max(calib['predicted'].max(), calib['observed'].max()) * 1.12
ax.plot([0, top], [0, top], color=INK2, linewidth=1, linestyle='--', label='Perfect calibration')
ax.plot(calib['predicted'], calib['observed'], color=BLUE, linewidth=2, marker='o', markersize=7,
        markeredgecolor=SURFACE, markeredgewidth=2, label='Final model')
ax.annotate('Final model', (calib['predicted'].iloc[-2], calib['observed'].iloc[-2]),
            xytext=(10, -16), textcoords='offset points', ha='left', color=INK)
ax.set_xlim(0, top)
ax.set_ylim(0, top)
ax.set_xlabel('Predicted risk (mean per tenth)')
ax.set_ylabel('Observed seropositive rate')
ax.set_title('Calibration on the test set (n = 1,002)', loc='left')
ax.legend(frameon=False, loc='upper left', labelcolor=INK2)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'calibration_final.png'), dpi=200)
plt.close(fig)


# =====================================================================
# STEP 7: Cut-off metrics and decision curve
# =====================================================================
flag = p_test >= cutoff
tp, fp = int((flag & (yt == 1)).sum()), int((flag & (yt == 0)).sum())
fn, tn = int((~flag & (yt == 1)).sum()), int((~flag & (yt == 0)).sum())
print(f'\nAt the training cut-off ({cutoff:.3f}) on the test set')
print(f'  Flagged {flag.sum()} of {len(yt)} people ({flag.mean():.0%})   TP {tp}  FP {fp}  FN {fn}  TN {tn}')
print(f'  Sensitivity {tp / (tp + fn):.3f}   Specificity {tn / (tn + fp):.3f}   '
      f'PPV {tp / (tp + fp):.3f}   NPV {tn / (tn + fn):.3f}')

# Net benefit = true positives found minus false positives, weighted by the risk threshold
pts = np.arange(0.05, 0.41, 0.01)
prev = yt.mean()
nb_model = [((p_test >= t) & (yt == 1)).mean() - ((p_test >= t) & (yt == 0)).mean() * t / (1 - t) for t in pts]
nb_all = [prev - (1 - prev) * t / (1 - t) for t in pts]
print('  Decision curve (net benefit per person):')
for t in [0.10, 0.15, 0.20, 0.25, 0.30]:
    i = int(round((t - 0.05) / 0.01))
    print(f'    threshold {t:.2f}: model {nb_model[i]:+.4f}   test everyone {nb_all[i]:+.4f}   test no one +0.0000')

fig, ax = plt.subplots(figsize=(5.6, 4.2))
ax.axhline(0, color=INK2, linewidth=1, label='Test no one')
ax.plot(pts, nb_all, color=ORANGE, linewidth=2, label='Test everyone')
ax.plot(pts, nb_model, color=BLUE, linewidth=2, label='Final model')
ax.annotate('Final model', (pts[12], nb_model[12]), xytext=(6, 6), textcoords='offset points', color=INK)
ax.annotate('Test everyone', (pts[8], nb_all[8]), xytext=(8, -6), textcoords='offset points', color=INK)
ax.set_ylim(-0.02, prev + 0.01)
ax.set_xlim(pts[0], pts[-1])
ax.set_xlabel('Risk threshold for testing')
ax.set_ylabel('Net benefit')
ax.set_title('Decision curve on the test set', loc='left')
ax.legend(frameon=False, loc='upper right', labelcolor=INK2)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'decision_curve_final.png'), dpi=200)
plt.close(fig)


# =====================================================================
# STEP 8 (RQ4): Adjusted odds ratios, training set, all Set C variables
# =====================================================================
# Unpenalised logistic regression fitted by Newton's method; 95% CIs from the standard errors.
d = train.copy()
for col in ['INDFMPIR', 'HSD010', 'drinks_per_day', 'phq9_total_score']:
    d[col + '_missing'] = d[col].isna().astype(float)          # "was missing" flags, not reported
    d[col] = d[col].fillna(d[col].median())
for col in ['DMDEDUC2', 'SLD010H']:
    d[col] = d[col].fillna(d[col].median())
for col in ['yrs_in_us_grp', 'SLQ050', 'smoking_status']:
    d[col] = d[col].fillna(d[col].mode()[0])

design = pd.DataFrame(index=d.index)
design['Age (per 10 years)'] = d['RIDAGEYR'] / 10
design['Male (vs female)'] = d['RIAGENDR'] == 1
for code, name in {1: 'Mexican American', 2: 'Other Hispanic', 4: 'Non-Hispanic Black',
                   6: 'Non-Hispanic Asian', 7: 'Other / multiracial'}.items():
    design[name + ' (vs NH White)'] = d['RIDRETH3'] == code
for code, name in {1: 'Foreign-born, <10 y in US', 2: 'Foreign-born, 10-19 y', 3: 'Foreign-born, 20+ y',
                   4: 'Foreign-born, years unknown'}.items():
    design[name + ' (vs US-born)'] = d['yrs_in_us_grp'] == code
design['Education (per level, 1-5)'] = d['DMDEDUC2']
design['Income-to-poverty ratio (per unit)'] = d['INDFMPIR']
design['Toxocara positive'] = d['toxocara_pos']
design['Self-rated health (per level, 5 = poor)'] = d['HSD010']
design['Sleep hours (per hour)'] = d['SLD010H']
design['Told doctor about trouble sleeping'] = d['SLQ050'] == 1
design['Former smoker (vs never)'] = d['smoking_status'] == 1
design['Current smoker (vs never)'] = d['smoking_status'] == 2
design['Drinks per drinking day (per drink)'] = d['drinks_per_day']
design['PHQ-9 score (per 5 points)'] = d['phq9_total_score'] / 5
reported = list(design.columns)
for col in ['INDFMPIR', 'HSD010', 'drinks_per_day', 'phq9_total_score']:
    design[col + '_missing'] = d[col + '_missing']

Xm = np.column_stack([np.ones(len(design)), design.to_numpy(dtype=float)])
yv = y_train.to_numpy(dtype=float)
beta = np.zeros(Xm.shape[1])
for _ in range(100):
    p = 1 / (1 + np.exp(-Xm @ beta))
    step = np.linalg.solve(Xm.T @ (Xm * (p * (1 - p))[:, None]), Xm.T @ (yv - p))
    beta += step
    if np.abs(step).max() < 1e-10:
        break
p = 1 / (1 + np.exp(-Xm @ beta))
se = np.sqrt(np.diag(np.linalg.inv(Xm.T @ (Xm * (p * (1 - p))[:, None]))))

print(f'\nAdjusted odds ratios (training set, n = {len(design)}; all variables in one model)')
print(f'  {"Variable":44} {"OR":>5}  {"95% CI":>13}  {"p":>7}')
for j, name in enumerate(reported, start=1):
    lo, hi = np.exp(beta[j] - 1.96 * se[j]), np.exp(beta[j] + 1.96 * se[j])
    pval = 2 * stats.norm.sf(abs(beta[j] / se[j]))
    print(f'  {name:44} {np.exp(beta[j]):5.2f}  {lo:5.2f} - {hi:5.2f}  {pval:7.4f}')


# =====================================================================
# STEP 9 (RQ4): SHAP values for the final model on the test set
# =====================================================================
# For logistic regression SHAP is exact: coefficient x (value - average value) on the log-odds scale.
prepare, model = final.named_steps['prepare'], final.named_steps['model']
Xtr_t, Xte_t = prepare.transform(train[FINAL]), prepare.transform(test[FINAL])
sv = shap.LinearExplainer(model, shap.maskers.Independent(Xtr_t, max_samples=1000))(Xte_t).values

# Add the one-hot and "was missing" columns back up to their original feature
names = [n.split('__', 1)[1] for n in prepare.get_feature_names_out()]
owner = []
for n in names:
    n = n.replace('missingindicator_', '')
    owner.append(next(f for f in FINAL if n == f or n.startswith(f + '_')))
grouped = pd.DataFrame(sv, columns=names).T.groupby(owner).sum().T
importance = grouped.abs().mean().sort_values()
print('\nMean |SHAP| on the test set (log-odds), by feature')
for f, v in importance[::-1].items():
    print(f'  {LABELS[f]:26} {v:.3f}')

fig, ax = plt.subplots(figsize=(5.6, 3.4))
labels = [LABELS[f] for f in importance.index]
ax.barh(labels, importance.values, color=BLUE, height=0.6, edgecolor=SURFACE, linewidth=2)
for i, v in enumerate(importance.values):
    ax.text(v + importance.max() * 0.02, i, f'{v:.2f}', va='center', color=INK2)
ax.set_xlim(0, importance.max() * 1.18)
ax.set_axisbelow(True)
ax.grid(axis='y', visible=False)
ax.set_xlabel('Mean |SHAP value| (log-odds)')
ax.set_title('What drives the final model (test set)', loc='left')
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'shap_final.png'), dpi=200)
plt.close(fig)

print(f'\nFigures saved in {FIG_DIR}: calibration_final.png, decision_curve_final.png, shap_final.png')
