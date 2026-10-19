"""
NHANES 2013-2014 -> adult Toxoplasma gondii IgG dataset (simplified)

Load 8 files -> merge on SEQN -> keep adults with a clear lab result ->
clean SAS zeros + Refused/Don't-know codes -> build smoking, alcohol,
education, years-in-US, PHQ-9 and Toxocara features -> create target -> drop raw/leaky
columns -> save. Survey design columns (WTMEC2YR, SDMVPSU, SDMVSTRA) and SEQN
are kept for inference/joins and must NOT be used as model features.
"""
import os
import numpy as np
import pandas as pd

RAW_DIR = os.path.join('data', 'raw', 'NHANES')
OUT_CSV = os.path.join('data', 'processed', 'nhanes_adult_toxo_clean.csv')
WRITE_CSV = True


def read(name, cols):
    return pd.read_sas(os.path.join(RAW_DIR, name), format='xport')[cols]


# ---- 1. Load and merge (inner join on SEQN = exact intersection) -------------
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
print(f"Participants present in all {len(tables)} files: {len(df)}")

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
}
codes.update({item: [7, 9] for item in phq_items})
for col, bad in codes.items():
    df.loc[df[col].isin(bad), col] = np.nan

# ---- 5. Feature engineering --------------------------------------------------
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

# ~8% of people have no PHQ-9 at all, and most of them also lack HSD010 and
# drinks_per_day (they skipped that whole block). They are KEPT: dropping them
# under-represents foreign-born participants. Missing values are imputed
# inside the training pipeline.
no_phq = df['phq9_total_score'].isna()
print(f"Kept {int(no_phq.sum())} participants with missing PHQ-9 "
      f"({df.loc[no_phq, 'HSD010'].isna().mean():.0%} of them also lack HSD010)")

# ---- 6. Target + drop raw / leaky columns ------------------------------------
df['target_toxo_pos'] = (df['SSTOXG'] >= 33).astype(int)
drop_cols = ['SSTOXG', 'SMQ020', 'SMQ040', 'ALQ101', 'ALQ110', 'ALQ120Q',
             'ALQ130', 'DMDEDUC3', 'DMDYRSUS', 'SSTCAR'] + phq_items
df = df.drop(columns=drop_cols)

# ---- 7. Check and save -------------------------------------------------------
print(f"\nFinal shape: {df.shape}")
print("\nClass balance:")
print(df['target_toxo_pos'].value_counts().sort_index().to_string())
print("\nColumns still containing missing values:")
print(df.isna().sum()[lambda s: s > 0].to_string())
print("\nValues in coded columns:")
for col in ['RIAGENDR', 'RIDRETH3', 'DMDBORN4', 'yrs_in_us_grp', 'DMDEDUC2', 'HSD010',
            'SLQ050', 'smoking_status', 'toxocara_pos']:
    print(f"  {col:15}: {sorted(df[col].dropna().unique().tolist())}")
print(df[['RIDAGEYR', 'INDFMPIR', 'SLD010H', 'drinks_per_day', 'phq9_total_score']]
      .describe().loc[['min', 'mean', 'max']].round(2).to_string())

if WRITE_CSV:
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nSaved: {OUT_CSV}")