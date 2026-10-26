# NHANES Toxoplasma Project: Summary

**Status:** wrapped up on 26 Oct 2026. The analysis is finished, and the held-out test set has been used once.
**Next project:** self-distillation on the Makassar dataset (see section 11).

---

## 1. In one paragraph
We tested how well survey data from **NHANES 2013–2014** (5,006 U.S. adults, 750 = 15.0% positive for *Toxoplasma gondii* IgG) can predict infection. Four models (logistic regression, random forest, gradient boosting, XGBoost) were compared on four nested feature sets with repeated cross-validation and a corrected t-test. No feature group beyond demographics and no complex model gave a meaningful gain. The final model is therefore **logistic regression on 6 demographic variables**, which needs no blood test. On the untouched test set it reached **ROC-AUC 0.700 (95% CI 0.654–0.747)**, PR-AUC 0.316 (baseline 0.150), and its risks were well calibrated (slope 0.90). Birthplace and years in the U.S., age and race/ethnicity drive the predictions. The ceiling (~0.70–0.72) comes from the data: NHANES does not ask about cats, undercooked meat or soil.

---

## 2. Status at a glance

| Item | Status |
|---|---|
| Data cleaning (v1, v2) | Done |
| Model comparison, cross-validation (H1–H5) | Done |
| Final model + test-set evaluation (RQ4, RQ5) | Done, test set used **once** (do not re-tune on it) |
| Proposal (Claude Doc) with full Results section | Done |
| IEEE LaTeX paper | Draft done (`paper/toxoplasma_ieee.tex`) |
| Title, author affiliations, venue/year of Ichwan et al. ref | **Open** |
| Git commit of the new work | **Open** (see section 10) |

---

## 3. Data

**Source:** CDC NHANES 2013–2014 public files in `data/raw/NHANES/`
- Base (8): `SSTOXO_H` (Toxoplasma IgG, the label), `DEMO_H`, `HSQ_H`, `SLQ_H`, `DPQ_H` (PHQ-9), `SMQ_H`, `ALQ_H`, `SSTOCA_H` (Toxocara)
- Added for RQ2 (5): `SXQ_H` (sexual behaviour), `DUQ_H` (drug use), `OHXDEN_H` (dental exam), `OHQ_H` (oral health questionnaire), `HIV_H`
- Download URL pattern: `https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/<NAME>.xpt`

**Cohort and label**
- Adults aged 18+ with an IgG result: **n = 5,006**
- Positive = IgG at or above **33 IU/mL**. Borderline results (27 to 33) are coded **negative**, as the CDC instructs (11 people).
- Split: stratified 80/20, `random_state=42`, giving **4,004 training (600 positive) / 1,002 test (150 positive)**

**Processed files**

| File | Made by | Contents |
|---|---|---|
| `data/processed/nhanes_adult_toxo_clean.csv` (v1) | `scripts/cleaning2.py` | 19 columns: ID, 3 survey-design, target, 14 features |
| `data/processed/nhanes_adult_toxo_clean_v2.csv` (v2) | `scripts/cleaning3.py` | v1 + 15 features matching Ichwan et al. (34 columns). The same people in the same order, and the v1 columns are identical. |

Survey design columns (`WTMEC2YR`, `SDMVPSU`, `SDMVSTRA`) and `SEQN` are **never** model features.

---

## 4. Scripts and how to rerun

Run from the project root, in this order:

| Step | Script | What it does | Output |
|---|---|---|---|
| 1 | `scripts/cleaning3.py` | Builds v2 (v1 pipeline + Ichwan-style features) | `nhanes_adult_toxo_clean_v2.csv` |
| 2 | `scripts/train_simple.py` | Cross-validation of all models x feature sets A–D, hypothesis tests H1–H5 (training set only) | Printed table |
| 3 | `scripts/final_evaluation.py` | Tuning, cut-off, **one-time test-set evaluation**, calibration, decision curve, odds ratios, SHAP | `outputs/final_evaluation_results.txt`, `outputs/figures/*_final.png` |

Other scripts (older or helper versions, kept for history):
- `cleaning2.py`: builds v1 (still valid; v2 contains it)
- `cleaning.py`, `01_build_cleaned_dataset.py`: earlier cleaning versions
- `train_baseline.py`: earlier, longer training script (uses v1)
- `check_toxo_positives_negatives.py`: verifies positive/negative counts
- `fetch_sources.py`: downloads PubMed/openFDA/UniProt material (from an earlier LLM idea, now out of scope)

Packages: pandas, numpy, scikit-learn 1.7, xgboost 3.4, shap 0.52, matplotlib, scipy (`statsmodels` is **not** needed).

---

## 5. Key methods and decisions (and why)

| Decision | Why |
|---|---|
| Borderline IgG (27–33) coded negative | CDC documentation: retested and should be coded negative |
| Keep people with missing questionnaire data; impute inside the pipeline + "was missing" flags | Dropping them removes many foreign-born adults |
| Skip logic: a question skipped because it did not apply = 0 (e.g. no partner in past year means condomless sex = 0) | Otherwise 400+ people wrongly missing |
| 5-fold CV x 3 repeats (15 folds), same folds for every model | Fair comparison without touching the test set |
| Corrected t-test (Nadeau–Bengio); "meaningful" = ROC-AUC gain >= 0.01 **and** p < 0.05 | Folds overlap, so a normal t-test is too optimistic |
| If a feature group or complex model fails that rule, **keep the simpler option** | Gives the final model: LR on demographics, no blood test |
| No SMOTE / class weights | Keeps real probabilities, so calibration can be checked |
| Cut-off chosen on out-of-fold **training** predictions (Youden) | The test set must not choose anything |
| Survey weights not used | It's a prediction model for this sample; stated as a limitation |

**Feature sets**
- **A** demographics: age, sex, race/ethnicity, years in U.S., education, income-to-poverty ratio
- **B** = A + Toxocara
- **C** = B + self-rated health, sleep hours, sleep trouble, smoking, drinks/day, PHQ-9
- **D** = C + 14 Ichwan-style features (sexual behaviour, drug use, oral health, HIV)

---

## 6. Results

### Cross-validated ROC-AUC (training set, mean of 15 folds)

| Model | A | B | C | D |
|---|---|---|---|---|
| Logistic regression | 0.720 | **0.724** | 0.719 | 0.714 |
| Random forest | 0.710 | 0.716 | 0.718 | 0.717 |
| Gradient boosting | 0.707 | 0.711 | 0.704 | 0.694 |
| XGBoost | 0.710 | 0.714 | 0.712 | 0.700 |

Baselines: dummy 0.500, age-only LR 0.632. Best CV: LR on set B (ROC 0.724, PR 0.337, Brier 0.116).

### Hypotheses
- **H1** models beat the baselines: **yes** (all, p < 0.001)
- **H2** trees beat LR: **no** (differences −0.016 to −0.001)
- **H3** Toxocara helps: consistent **+0.004 to +0.006**, but below the 0.01 rule
- **H4** behaviour/health helps: **no** (−0.007 to +0.003)
- **H5** Ichwan-style features help: **no** (−0.012 to −0.001; sexual features slightly hurt XGBoost)

### Final model on the test set (LR, set A, C = 1; tuning gave no gain)

| Model | ROC-AUC (95% CI) | PR-AUC (95% CI) | Brier |
|---|---|---|---|
| Dummy | 0.500 | 0.150 | 0.127 |
| Age-only LR | 0.641 (0.593–0.689) | 0.216 (0.175–0.266) | 0.124 |
| **Final model** | **0.700 (0.654–0.747)** | **0.316 (0.253–0.393)** | **0.118** |

- Gain over age alone: **+0.059 (95% CI +0.018 to +0.103)**
- Calibration: mean predicted 14.9% vs observed 15.0%; slope 0.90, intercept −0.15
- At cut-off 0.124 (chosen on training data): flags 47%, **sensitivity 0.73, specificity 0.57, PPV 0.23, NPV 0.92**
- Decision curve: the model beats "test everyone" and "test no one" for thresholds 10–40%, but the gain is small

### Adjusted odds ratios (training set, all set C variables)
Age per 10 y **1.34** · male **1.26** · Other Hispanic **1.55** · Mexican American **0.64** · NH Asian **0.34** · foreign-born <10 y **6.13**, 10–19 y **3.26**, 20+ y **2.54**, unknown **2.28** · education per level **0.85** · **Toxocara 1.91** (linked to infection but barely predictive, because only ~7% have it) · income, health, sleep, smoking, alcohol, PHQ-9: not significant.

### SHAP (what drives the final model)
Birthplace / years in U.S. 0.48 > age 0.45 > race/ethnicity 0.30 > education 0.18 > sex 0.13 > income 0.08

---

## 7. Where everything is

| What | Where |
|---|---|
| Full test-set output | `outputs/final_evaluation_results.txt` |
| Figures | `outputs/figures/calibration_final.png`, `decision_curve_final.png`, `shap_final.png` (copies also in `paper/`) |
| IEEE LaTeX paper (for Prism) | `paper/toxoplasma_ieee.tex` (upload the 3 PNGs with it) |
| Revised proposal (Claude Doc, 25 refs, full Results) | https://claude.ai/artifact/7afvqusm9Rpopb9tc6775H |
| Original 34-item reference list | `paper/references.txt` (includes the old LLM references) |

---

## 8. Papers reviewed (notes)

- **Ichwan et al.** (HIV cohort, Indonesia): 78 people after filtering, SMOTE + PCA, one 70/30 split (24 test), best weighted F1 0.80 in one subgroup. Venue/year still unknown. Our v2 tested whether its features carry over: they **don't** in the U.S. general population.
- **Guo et al. 2024** (vitamin D, NHANES): used as a method template only. Its AUC of 1.0 is not believable; its own LR/GBM reached 0.77–0.79.
- **Team draft `main.pdf`** (pooled 2009–2014, 44 biological markers, best test ROC-AUC 0.730). Stronger analysis and a good main message ("the limit is information, not the learner"). Concerns to fix:
  1. Check that the **"associated markers" selection happens inside each CV fold** (otherwise CV 0.740 is optimistic)
  2. Blood markers clash with the "survey answers" motivation; add a survey-only model
  3. Do not lean on test-set gains of 0.015–0.020 (within noise)
  4. Ref [15] (Polanunu 2021) does **not** report any ROC-AUC. Replace the "ROC-AUCs above 0.80" sentence (our Makassar analysis gives ~0.87 and can be cited as our own analysis of their data)
  5. Missing: ethics statement, recent references, the Ichwan-feature test, decision curve

---

## 9. Other datasets checked

| Dataset | Verdict |
|---|---|
| NHANES 2009–2010, 2011–2012 | **Same test and cut-off as 2013–2014, can be pooled** (6,670 and 6,097 with IgG). 2011–12 includes special Toxoplasma subsample weights (`WTSSTT2Y`). |
| NHANES 1999–2004 | Ages 6–49 only, different calibration (positive at 10 IU/mL or above), so not directly comparable |
| NHANES 2005–2008, 2015 onwards | No Toxoplasma testing |
| **Makassar (Polanunu et al. 2021, PLoS ONE)** | Public, 184 pregnant women, 60 positive, **direct exposure questions**. Saved at `data/raw/makassar/polanunu2021_S1.sav`. Exposures alone give CV ROC-AUC ~0.86–0.88 vs 0.50 for demographics. |
| UK Biobank | ~9,400 with Toxoplasma serology, but needs an approved application + fees |
| Flegr (`data/raw/flegr/`) | Only 747 labelled rows, 10 columns, self-reported |
| Achtziger (`data/raw/achtziger/`) | Clinical study: 13 positives. Online survey: no real label |
| `data/Data Clean.csv` | Indonesian health survey, no Toxoplasma variable |

---

## 10. Open items
- [ ] Decide the paper title (only demographics stay in the final model)
- [ ] Fill in author department, university and emails in the LaTeX
- [ ] Find the venue and year for Ichwan et al.
- [ ] Commit the new work to git. Untracked: `scripts/cleaning3.py`, `scripts/train_simple.py`, `scripts/final_evaluation.py`, `data/processed/nhanes_adult_toxo_clean_v2.csv`, the 5 new NHANES files, `data/raw/makassar/`, `paper/`, `outputs/`, this file
- [ ] (Team draft) Check the marker-selection leakage in the code

---

## 11. Next project: self-distillation on the Makassar data
Plan agreed so far (from Mr Feri's suggestion: KD, MTKD, deep mutual learning, self-distillation; self-distillation is the cheapest):
1. `scripts/clean_makassar.py` → `data/processed/makassar_clean.csv`
2. Baselines (LR, RF, XGBoost, small MLP) with repeated CV
3. Controls: label smoothing, stronger L2
4. Self-distillation with an **out-of-fold teacher**; LR/XGBoost learn soft labels via the duplicate-row + sample-weight trick, the MLP learns them directly
5. Optional: MTKD (ensemble teacher), deep mutual learning
6. Size experiment: compare with NHANES subsamples (184, 500, 1,000, 4,004)

**Expectations:** ROC-AUC likely unchanged (the signal is simple and strong); calibration and stability are where gains are possible; a null result is still publishable if framed well.
