"""
NHANES 2013-2014 -> adult Toxoplasma gondii IgG dataset, v2 (IEEE-paper features)

Same pipeline as cleaning2.py (same cohort, same base features), plus features
that mirror the questionnaire columns of the IEEE paper "Machine Learning-Based
Prediction of Seropositive Toxoplasma gondii from Lifestyle Factors":
  - sexual behaviour      (SXQ_H)   -> condom use, oral-sex protection, anal sex, partners, STDs
  - drug use              (DUQ_H)   -> injection drug use, hard drug use
  - oral health           (OHXDEN_H, OHQ_H) -> missing teeth, loose teeth, mouth aching
  - HIV antibody          (HIV_H)

The 8 base files are inner-joined (identical cohort to cleaning2.py). The new
files are LEFT-joined: they were only asked of narrower age ranges (SXQ/DUQ
18-69, many items 18-59; HIV 18-59; OHQ620/855 30+), so an inner join would
silently drop older adults. Not-asked = NaN, imputed inside the training pipeline.

Writes data/processed/nhanes_adult_toxo_clean_v2.csv (the v1 file is untouched).
Survey design columns (WTMEC2YR, SDMVPSU, SDMVSTRA) and SEQN are kept for
inference/joins and must NOT be used as model features.
"""
import os
import numpy as np
import pandas as pd

RAW_DIR = os.path.join('data', 'raw', 'NHANES')
OUT_CSV = os.path.join('data', 'processed', 'nhanes_adult_toxo_clean_v2.csv')
WRITE_CSV = True


def read(name, cols):
    return pd.read_sas(os.path.join(RAW_DIR, name), format='xport')[cols]


# ---- 1. Load and merge -------------------------------------------------------
phq_items = [f'DPQ0{i}0' for i in range(1, 10)]          # DPQ010 ... DPQ090
tables = [
    read('SSTOXO_H.xpt', ['SEQN', 'SSTOXG']),
    read('DEMO_H.xpt', ['SEQN', 'RIDAGEYR', 'RIAGENDR', 'RIDRETH3', 'DMDBORN4',
                        'DMDYRSUS', 'DMDEDUC2', 'DMDEDUC3', 'INDFMPIR',
                        'WTMEC2YR', 'SDMVPSU', 'SDMVSTRA']),     # survey design kept for inference
    read('HSQ_H.xpt', ['SEQN', 'HSD010']),
    read('SLQ_H.xpt', ['SEQN', 'SLD010H', 'SLQ050']),
    read('DPQ_H.xpt', ['SEQN'] + phq_items),
    read('SMQ_H.xpt', ['SEQN', 'SMQ020', 'SMQ040']),
    read('ALQ_H.xpt', ['SEQN', 'ALQ101', 'ALQ110', 'ALQ120Q', 'ALQ130']),
    read('SSTOCA_H.xpt', ['SEQN', 'SSTCAR']),     # Toxocara antibody: exposure marker (soil / dog-cat feces)
]
df = tables[0]
for t in tables[1:]:
    df = df.merge(t, on='SEQN', how='inner')
print(f"Participants present in all {len(tables)} base files: {len(df)}")

sxq_cols = ['SXD021', 'SXQ803', 'SXQ806', 'SXQ809', 'SXQ853',     # men
            'SXQ703', 'SXQ706', 'SXQ709', 'SXQ741',               # women
            'SXQ645', 'SXQ610', 'SXQ251', 'SXD031', 'SXQ836',
            'SXD171', 'SXQ410', 'SXD101', 'SXQ130',
            'SXD510', 'SXQ550', 'SXD450', 'SXQ490',               # past-year partners
            'SXQ260', 'SXQ265', 'SXQ270', 'SXQ272']
duq_cols = ['DUQ240', 'DUQ370', 'DUQ410']
tooth_cols = [f'OHX{i:02d}TC' for i in range(1, 33)]
ohq_cols = ['OHQ620', 'OHQ855']
extra = [
    read('SXQ_H.xpt', ['SEQN'] + sxq_cols),
    read('DUQ_H.xpt', ['SEQN'] + duq_cols),
    read('OHXDEN_H.xpt', ['SEQN', 'OHDDESTS'] + tooth_cols),
    read('OHQ_H.xpt', ['SEQN'] + ohq_cols),
    read('HIV_H.xpt', ['SEQN', 'LBDHI']),
]
for t in extra:
    df = df.merge(t, on='SEQN', how='left')

# ---- 2. Cohort: adults with a clear lab result -------------------------------
df = df[(df['RIDAGEYR'] >= 18) & df['SSTOXG'].notna()].copy()
# CDC SSTOXO_H documentation: "All results originally considered Equivocal were
# repeated at least twice for confirmation and should be coded as Negative."
# So 27 <= SSTOXG < 33 is kept and falls on the negative side of the >= 33 cutoff.
equivocal = (df['SSTOXG'] >= 27) & (df['SSTOXG'] < 33)
print(f"Kept {int(equivocal.sum())} equivocal samples (27 <= SSTOXG < 33) as Negative, per CDC")
print(f"Adult cohort: {len(df)}")

# ---- 3. SAS float artifact: pandas reads a stored 0 as 5.3976e-79 ------------
for col in df.select_dtypes('float').columns:
    tiny = (df[col].abs() > 0) & (df[col].abs() < 1e-70)
    df.loc[tiny, col] = 0.0

# ---- 4. Refused / Don't-know codes -> NaN ------------------------------------
# (SLD010H: a plain 7 is a valid "7 hours", so only 77/99 are removed)
codes = {
    'DMDBORN4': [77, 99], 'DMDYRSUS': [77, 99], 'DMDEDUC2': [7, 9], 'DMDEDUC3': [77, 99],
    'HSD010': [7, 9], 'SLD010H': [77, 99], 'SLQ050': [7, 9],
    'SMQ020': [7, 9], 'SMQ040': [7, 9],
    'ALQ101': [7, 9], 'ALQ110': [7, 9], 'ALQ120Q': [777, 999],
    'ALQ130': [77, 99, 777, 999],
    # SXQ645: 5 = "Unsure" is not a frequency either
    'SXQ645': [5, 7, 9], 'SXQ610': [77, 99], 'SXD031': [77, 99],
    'DUQ410': [77, 99], 'LBDHI': [3],                              # 3 = indeterminate
}
codes.update({item: [7, 9] for item in phq_items})
codes.update({c: [7, 9] for c in ['SXD021', 'SXQ803', 'SXQ806', 'SXQ809', 'SXQ853',
                                  'SXQ703', 'SXQ706', 'SXQ709', 'SXQ741', 'SXQ251',
                                  'SXQ260', 'SXQ265', 'SXQ270', 'SXQ272',
                                  'DUQ240', 'DUQ370', 'OHQ620', 'OHQ855']})
codes.update({c: [77777, 99999] for c in ['SXQ836', 'SXD171', 'SXQ410', 'SXD101', 'SXQ130',
                                          'SXD510', 'SXQ550', 'SXD450', 'SXQ490']})
for col, bad in codes.items():
    df.loc[df[col].isin(bad), col] = np.nan

# ---- 5. Feature engineering (base, as in cleaning2.py) -----------------------
# Smoking: 0 = never, 1 = former, 2 = current
df['smoking_status'] = np.nan
df.loc[df['SMQ020'] == 2, 'smoking_status'] = 0
df.loc[(df['SMQ020'] == 1) & (df['SMQ040'] == 3), 'smoking_status'] = 1
df.loc[(df['SMQ020'] == 1) & df['SMQ040'].isin([1, 2]), 'smoking_status'] = 2

# Alcohol: drinks per drinking day; 0 for lifetime abstainers and anyone
# who did not drink at all in the past 12 months
df['drinks_per_day'] = df['ALQ130']
non_drinker = ((df['ALQ101'] == 2) & (df['ALQ110'] == 2)) | (df['ALQ120Q'] == 0)
df.loc[non_drinker, 'drinks_per_day'] = 0.0

# Years in US (DMDYRSUS is only asked of foreign-born):
# 0 = US-born, 1 = <10 yrs, 2 = 10-19 yrs, 3 = 20+ yrs, 4 = foreign-born, unknown
df['yrs_in_us_grp'] = np.nan
df.loc[df['DMDBORN4'] == 1, 'yrs_in_us_grp'] = 0
df.loc[df['DMDYRSUS'].isin([1, 2, 3]), 'yrs_in_us_grp'] = 1
df.loc[df['DMDYRSUS'].isin([4, 5]), 'yrs_in_us_grp'] = 2
df.loc[df['DMDYRSUS'].isin([6, 7, 8, 9]), 'yrs_in_us_grp'] = 3
df.loc[(df['DMDBORN4'] == 2) & df['DMDYRSUS'].isna(), 'yrs_in_us_grp'] = 4

# Toxocara antibody (SSTCAR: 1 = Positive, 2 = Negative) -> 1/0
df['toxocara_pos'] = df['SSTCAR'].map({1.0: 1, 2.0: 0})


# Education: 18-19 year olds answer DMDEDUC3 instead of DMDEDUC2
def edu3_to_edu2(x):
    if pd.isna(x):
        return np.nan
    if x <= 8 or x in (55, 66):      # less than 9th grade
        return 1.0
    if x <= 12:                      # 9th-12th grade, no diploma
        return 2.0
    if x in (13, 14):                # high school graduate / GED
        return 3.0
    if x == 15:                      # more than high school
        return 4.0
    return np.nan


young = df['DMDEDUC2'].isna() & df['RIDAGEYR'].isin([18, 19])
df.loc[young, 'DMDEDUC2'] = df.loc[young, 'DMDEDUC3'].map(edu3_to_edu2)

# PHQ-9 total (0-27); NaN if any item is missing
df['phq9_total_score'] = df[phq_items].sum(axis=1, skipna=False)

no_phq = df['phq9_total_score'].isna()
print(f"Kept {int(no_phq.sum())} participants with missing PHQ-9 "
      f"({df.loc[no_phq, 'HSD010'].isna().mean():.0%} of them also lack HSD010)")

# ---- 6. Feature engineering (IEEE-paper features) ----------------------------
male = df['RIAGENDR'] == 1
female = df['RIAGENDR'] == 2
never_sex = df['SXD021'] == 2       # skipped the rest of SXQ; their "exposure" is 0


def yes_no(s):
    return s.map({1.0: 1, 2.0: 0})


# Paper: "Frequency of anal intercourse" -> NHANES only asks "ever"
# Men: anal sex with a woman (SXQ806) or >= 1 male anal partner (SXQ836)
# Women: anal sex with a man (SXQ706)
anal = pd.concat([yes_no(df['SXQ806']).where(male),
                  (df['SXQ836'] > 0).astype(float).where(male & df['SXQ836'].notna()),
                  yes_no(df['SXQ706']).where(female)], axis=1)
df['anal_sex_ever'] = anal.max(axis=1)          # 1 if any source says yes, NaN if none answered
df.loc[never_sex, 'anal_sex_ever'] = 0

# Same-sex partner ever (men: SXQ809, women: SXQ709); the paper's cohort was MSM
df['same_sex_partner_ever'] = yes_no(df['SXQ809']).where(male).fillna(yes_no(df['SXQ709']).where(female))
df.loc[never_sex, 'same_sex_partner_ever'] = 0

# Lifetime partners: opposite-sex count + same-sex count (same-sex question is
# only asked of people who reported a same-sex partner, so not-asked = 0)
opp = df['SXD171'].where(male).fillna(df['SXD101'].where(female))
same = df['SXQ410'].where(male).fillna(df['SXQ130'].where(female))
same = same.where(df['same_sex_partner_ever'] != 0, 0)
df['lifetime_partners'] = opp + same.fillna(0)
df.loc[never_sex, 'lifetime_partners'] = 0

# Vaginal/anal sex in past 12 months (SXQ610): 0 = never ... 6 = 365+ times.
# SXQ610/SXQ251 are skipped for people with no partner in the past year -> 0
opp_year = df['SXD510'].where(male).fillna(df['SXD450'].where(female))
same_year = df['SXQ550'].where(male).fillna(df['SXQ490'].where(female))
same_year = same_year.where(df['same_sex_partner_ever'] != 0, 0)
no_partner_year = (opp_year + same_year.fillna(0)) == 0
df['sex_freq_year'] = df['SXQ610']
df.loc[never_sex | no_partner_year, 'sex_freq_year'] = 0

# Paper: "Partner's condom use" (higher = less use) -> SXQ251 "how often sex
# WITHOUT a condom": 1 = never ... 5 = always  ->  0 ... 4
df['condomless_sex_freq'] = df['SXQ251'] - 1
df.loc[never_sex | (df['sex_freq_year'] == 0), 'condomless_sex_freq'] = 0

# Paper: "Condom use during oral sex" (higher = less use) -> SXQ645
# "use protection when performing oral sex": 1 = never ... 4 = always  ->  3 ... 0
df['oral_sex_no_protection'] = 4 - df['SXQ645']
gave_oral = pd.concat([df['SXQ803'].where(male), df['SXQ853'].where(male),
                       df['SXQ703'].where(female), df['SXQ741'].where(female)], axis=1)
never_gave_oral = (gave_oral == 2).any(axis=1) & ~(gave_oral == 1).any(axis=1)
df.loc[never_sex | never_gave_oral, 'oral_sex_no_protection'] = 0

df['age_first_sex'] = df['SXD031']                          # 9 = "9 or younger"

# Any doctor-diagnosed STD: herpes, genital warts, gonorrhea, chlamydia
std = df[['SXQ260', 'SXQ265', 'SXQ270', 'SXQ272']]
df['std_ever'] = np.where((std == 1).any(axis=1), 1.0, np.where((std == 2).all(axis=1), 0.0, np.nan))

# Paper: "Current injection drug use" / "History of used needle use" (DUQ_H)
df['inject_drugs_ever'] = yes_no(df['DUQ370'])
# Lifetime injection count category: 0 = never ... 6 = 100+ times
df['inject_count_grp'] = df['DUQ410']
df.loc[df['DUQ370'] == 2, 'inject_count_grp'] = 0
df['hard_drug_ever'] = yes_no(df['DUQ240'])                 # cocaine / heroin / meth

# Paper: "Tooth loss" -> missing permanent teeth from the dental exam, 28-tooth
# count (third molars #1, 16, 17, 32 excluded: often never erupt or are removed
# electively). Missing = 4 "not present" or 3 "implant". Complete exams only.
teeth28 = [f'OHX{i:02d}TC' for i in range(1, 33) if i not in (1, 16, 17, 32)]
df['missing_teeth'] = df[teeth28].isin([3, 4]).sum(axis=1).astype(float)
df.loc[df['OHDDESTS'] != 1, 'missing_teeth'] = np.nan
# Paper: "Tooth loss due to natural causes" -> teeth loosened without injury (30+)
df['loose_teeth'] = yes_no(df['OHQ855'])
# Paper: "Mouth ulcers" (closest available) -> aching in mouth past year (30+)
# OHQ620: 1 = very often ... 5 = never  ->  4 ... 0
df['mouth_aching_freq'] = 5 - df['OHQ620']

df['hiv_pos'] = df['LBDHI'].map({1.0: 1, 2.0: 0})

# ---- 7. Target + drop raw / leaky columns ------------------------------------
df['target_toxo_pos'] = (df['SSTOXG'] >= 33).astype(int)
drop_cols = (['SSTOXG', 'SMQ020', 'SMQ040', 'ALQ101', 'ALQ110', 'ALQ120Q',
              'ALQ130', 'DMDEDUC3', 'DMDYRSUS', 'SSTCAR', 'OHDDESTS', 'LBDHI']
             + phq_items + sxq_cols + duq_cols + tooth_cols + ohq_cols)
df = df.drop(columns=drop_cols)

# ---- 8. Check and save -------------------------------------------------------
new_cols = ['anal_sex_ever', 'same_sex_partner_ever', 'lifetime_partners', 'sex_freq_year',
            'condomless_sex_freq', 'oral_sex_no_protection', 'age_first_sex', 'std_ever',
            'inject_drugs_ever', 'inject_count_grp', 'hard_drug_ever',
            'missing_teeth', 'loose_teeth', 'mouth_aching_freq', 'hiv_pos']
print(f"\nFinal shape: {df.shape}")
print("\nClass balance:")
print(df['target_toxo_pos'].value_counts().sort_index().to_string())
print("\nColumns still containing missing values:")
print(df.isna().sum()[lambda s: s > 0].to_string())
print("\nValues in coded columns:")
for col in ['RIAGENDR', 'RIDRETH3', 'DMDBORN4', 'yrs_in_us_grp', 'DMDEDUC2', 'HSD010',
            'SLQ050', 'smoking_status', 'toxocara_pos'] + new_cols:
    vals = sorted(df[col].dropna().unique().tolist())
    print(f"  {col:22}: {vals if len(vals) <= 10 else f'{len(vals)} values, {vals[0]} .. {vals[-1]}'}")
print(df[['RIDAGEYR', 'INDFMPIR', 'SLD010H', 'drinks_per_day', 'phq9_total_score',
          'lifetime_partners', 'missing_teeth']]
      .describe().loc[['min', 'mean', 'max']].round(2).to_string())

print("\nShare missing in new features, by age band:")
age_band = pd.cut(df['RIDAGEYR'], [17, 29, 59, 69, 120], labels=['18-29', '30-59', '60-69', '70+'])
with pd.option_context('display.width', 200, 'display.max_columns', None):
    print(df[new_cols].isna().groupby(age_band, observed=True).mean().round(2).T.to_string())

if WRITE_CSV:
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nSaved: {OUT_CSV}")
