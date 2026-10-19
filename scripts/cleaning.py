"""
NHANES 2013-2014 -> Adult Toxoplasma gondii IgG seropositivity dataset
================================================================================
WORKFLOW:
 1. Load 7 official CDC NHANES .xpt files from data/raw/NHANES/
 2. Find the exact participant intersection on 'SEQN'
 3. Define the cohort:
      - adults only (RIDAGEYR >= 18)
      - drop participants with no SSTOXG result (never label NaN as negative)
      - drop equivocal titers (27 <= SSTOXG < 33 IU/mL)
 4. Clean SAS float artifacts (5.3976e-79 -> 0.0) in ALL float columns
 5. Clean special CDC codes (Refused / Don't know) -> NaN
 6. Resolve survey skip logic:
      - Smoking   (SMQ020 + SMQ040                  -> smoking_status 0/1/2)
      - Alcohol   (ALQ101 + ALQ110 + ALQ120Q + ALQ130 -> drinks_per_day)
      - Education for 18-19 year olds (DMDEDUC3 -> DMDEDUC2 scale)
 7. Compute composite clinical score:
      - phq9_total_score (sum of 9 depression items, 0 to 27)
 8. Define binary target:
      - target_toxo_pos = 1 if SSTOXG >= 33.0 IU/mL, else 0
 9. Drop raw lab titer (SSTOXG) and raw helper columns to prevent leakage
10. Verify, and (optionally) save to data/processed/nhanes_adult_toxo_clean.csv
    -> CSV writing is controlled by WRITE_CSV below (off until you've
       verified the printed results).
================================================================================
"""

import os
from functools import reduce

import numpy as np
import pandas as pd

# ------------------------------------------------------------------ settings
RAW_DIR = os.path.join('data', 'raw', 'NHANES')
PROCESSED_DIR = os.path.join('data', 'processed')
OUTPUT_CSV = os.path.join(PROCESSED_DIR, 'nhanes_adult_toxo_clean.csv')

WRITE_CSV = False        # flip to True once you're happy with the printed output

MIN_AGE = 18
NEG_CUTOFF = 27.0        # SSTOXG <  27         -> Negative (CDC)
POS_CUTOFF = 33.0        # SSTOXG >= 33         -> Positive (CDC)
                         # 27 <= SSTOXG < 33    -> Equivocal (excluded)
SAS_ZERO_TOL = 1e-70     # |x| below this (but > 0) is the XPT "zero" artifact


# ------------------------------------------------------------------- helpers
def banner(title):
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def read_xpt(filename, cols):
    """Read a SAS XPORT file and keep only the requested columns."""
    print(f"-> Reading {filename} ...")
    return pd.read_sas(os.path.join(RAW_DIR, filename), format='xport')[cols]


def map_edu3_to_edu2(x):
    """
    Map DMDEDUC3 (education, ages 6-19) onto the DMDEDUC2 (adults 20+) scale:
        1 = less than 9th grade
        2 = 9th-11th grade (includes 12th grade, no diploma)
        3 = high school graduate / GED
        4 = some college or AA degree
        5 = college graduate or above
    DMDEDUC3 codes: 0-12 = grade, 13 = HS grad, 14 = GED, 15 = more than HS,
                    55 = <5th grade, 66 = <9th grade, 77/99 = refused/unknown.
    """
    if pd.isna(x) or x in (77, 99):
        return np.nan
    if x <= 8 or x in (55, 66):
        return 1.0
    if x <= 12:
        return 2.0
    if x in (13, 14):
        return 3.0
    if x == 15:
        return 4.0
    return np.nan


# ---------------------------------------------------------------------- main
def main():
    # ------------------------------------------------------------------ 1
    banner("STEP 1: LOADING RAW CDC NHANES FILES & EXTRACTING COLUMNS")

    df_toxo = read_xpt('SSTOXO_H.xpt', ['SEQN', 'SSTOXG'])

    demo_cols = ['SEQN', 'RIDAGEYR', 'RIAGENDR', 'RIDRETH3', 'DMDBORN4',
                 'DMDEDUC2', 'DMDEDUC3', 'INDFMPIR']
    df_demo = read_xpt('DEMO_H.xpt', demo_cols)

    df_hsq = read_xpt('HSQ_H.xpt', ['SEQN', 'HSD010'])
    df_slq = read_xpt('SLQ_H.xpt', ['SEQN', 'SLD010H', 'SLQ050'])

    dpq_items = [f'DPQ0{i}0' for i in range(1, 10)]   # DPQ010 ... DPQ090
    df_dpq = read_xpt('DPQ_H.xpt', ['SEQN'] + dpq_items)

    df_smq = read_xpt('SMQ_H.xpt', ['SEQN', 'SMQ020', 'SMQ040'])
    df_alq = read_xpt('ALQ_H.xpt', ['SEQN', 'ALQ101', 'ALQ110', 'ALQ120Q', 'ALQ130'])

    # ------------------------------------------------------------------ 2
    banner("STEP 2: FINDING THE COMMON INTERSECTION ON 'SEQN'")

    frames = [df_toxo, df_demo, df_hsq, df_slq, df_dpq, df_smq, df_alq]
    common_seqn = set.intersection(*[set(f['SEQN']) for f in frames])
    print(f"Participants present in all {len(frames)} files: {len(common_seqn)}")

    merged = reduce(lambda left, right: left.merge(right, on='SEQN', how='inner'), frames)
    print(f"Merged shape: {merged.shape} (rows: {merged.shape[0]}, cols: {merged.shape[1]})")

    # ------------------------------------------------------------------ 3
    banner("STEP 3: DEFINING THE COHORT (ADULTS, VALID & NON-EQUIVOCAL SSTOXG)")

    n_start = len(merged)
    merged = merged[merged['RIDAGEYR'] >= MIN_AGE]
    print(f"  Age >= {MIN_AGE}: kept {len(merged)} of {n_start}")

    n_before = len(merged)
    merged = merged[merged['SSTOXG'].notna()]
    print(f"  Dropped {n_before - len(merged)} participants with missing SSTOXG")

    equivocal = (merged['SSTOXG'] >= NEG_CUTOFF) & (merged['SSTOXG'] < POS_CUTOFF)
    print(f"  Dropped {int(equivocal.sum())} equivocal samples "
          f"({NEG_CUTOFF} <= SSTOXG < {POS_CUTOFF} IU/mL)")
    merged = merged[~equivocal].copy()

    print(f"  Final cohort size: {len(merged)}")
    print(f"  Age range in cohort: {merged['RIDAGEYR'].min():.0f} - "
          f"{merged['RIDAGEYR'].max():.0f}")

    # ------------------------------------------------------------------ 4
    banner("STEP 4: CLEANING SAS FLOAT ZEROES (5.3976e-79 -> 0.0), ALL FLOAT COLUMNS")

    # pandas decodes a stored XPT zero as 5.397605346934028e-79. This matters for
    # columns where 0 is a real answer (PHQ-9 items, INDFMPIR, ALQ120Q, DMDEDUC3 ...),
    # so we fix every float column instead of a hand-picked list.
    for col in merged.select_dtypes(include='float').columns:
        abs_vals = merged[col].abs()
        mask = (abs_vals > 0) & (abs_vals < SAS_ZERO_TOL)
        n_fixed = int(mask.sum())
        if n_fixed > 0:
            merged.loc[mask, col] = 0.0
            print(f"  Fixed {n_fixed:4d} SAS zeroes in column: {col}")

    # ------------------------------------------------------------------ 5
    banner("STEP 5: CLEANING CDC SPECIAL CODES (Refused / Don't know -> NaN)")

    # Most questionnaire items use 7 = Refused, 9 = Don't know.
    # Some variables use 77/99 or 777/999 instead; hours of sleep (SLD010H) must
    # NOT lose its legitimate "7" (= 7 hours).
    special_codes = {
        'DMDBORN4': [77, 99],
        'DMDEDUC2': [7, 9],
        'DMDEDUC3': [77, 99],
        'HSD010':   [7, 9],
        'SLD010H':  [77, 99],
        'SLQ050':   [7, 9],
        'SMQ020':   [7, 9],
        'SMQ040':   [7, 9],
        'ALQ101':   [7, 9],
        'ALQ110':   [7, 9],
        'ALQ120Q':  [777, 999],
        'ALQ130':   [77, 99, 777, 999],
    }
    for item in dpq_items:
        special_codes[item] = [7, 9]

    for col, codes in special_codes.items():
        mask = merged[col].isin(codes)
        n_found = int(mask.sum())
        if n_found > 0:
            merged.loc[mask, col] = np.nan
            print(f"  Converted {n_found:4d} Refused/Don't-know entries to NaN in: {col}")

    # ------------------------------------------------------------------ 6
    banner("STEP 6: RESOLVING SURVEY SKIP LOGIC & FEATURE ENGINEERING")

    # A. SMOKING STATUS --------------------------------------------------
    # SMQ020: 1 = smoked >= 100 cigarettes in life, 2 = no
    # SMQ040: 1 = every day, 2 = some days, 3 = not at all (asked only if SMQ020 == 1)
    print("-> Engineering 'smoking_status' (0=Never, 1=Former, 2=Current)...")
    smoking_status = pd.Series(np.nan, index=merged.index)
    smoking_status[merged['SMQ020'] == 2] = 0.0
    smoking_status[(merged['SMQ020'] == 1) & (merged['SMQ040'] == 3)] = 1.0
    smoking_status[(merged['SMQ020'] == 1) & (merged['SMQ040'].isin([1, 2]))] = 2.0
    merged['smoking_status'] = smoking_status

    # B. ALCOHOL DRINKS PER DAY ------------------------------------------
    # Skip pattern (NHANES 2013-14):
    #   ALQ101 = 1                -> ALQ120Q (days drinking, past 12 months)
    #   ALQ101 = 2                -> ALQ110 (>= 12 drinks in lifetime?)
    #       ALQ110 = 1            -> ALQ120Q
    #       ALQ110 = 2 / 7 / 9    -> end of section
    #   ALQ120Q = 0               -> did not drink in the past year (ALQ130 skipped)
    #   ALQ130 = average drinks per drinking day over the past 12 months
    print("-> Engineering 'drinks_per_day'...")
    drinks_per_day = merged['ALQ130'].copy()
    lifetime_abstainer = (merged['ALQ101'] == 2) & (merged['ALQ110'] == 2)
    no_drinking_past_year = merged['ALQ120Q'] == 0
    zero_mask = lifetime_abstainer | no_drinking_past_year
    drinks_per_day[zero_mask] = 0.0
    merged['drinks_per_day'] = drinks_per_day
    print(f"  Set {int(zero_mask.sum())} non-drinkers (lifetime or past year) to 0.0")

    # C. EDUCATION FOR 18-19 YEAR OLDS -----------------------------------
    # DMDEDUC2 is only asked for ages 20+; 18-19 year olds answer DMDEDUC3.
    print("-> Resolving education for 18-19 year olds via DMDEDUC3...")
    needs_edu = merged['DMDEDUC2'].isna() & merged['RIDAGEYR'].isin([18, 19])
    filled = merged.loc[needs_edu, 'DMDEDUC3'].map(map_edu3_to_edu2)
    merged.loc[needs_edu, 'DMDEDUC2'] = filled
    print(f"  18-19 year olds needing education: {int(needs_edu.sum())}, "
          f"successfully mapped: {int(filled.notna().sum())}")

    # D. PHQ-9 TOTAL SCORE -----------------------------------------------
    print("-> Engineering 'phq9_total_score' (sum of 9 items, 0 to 27)...")
    merged['phq9_total_score'] = merged[dpq_items].sum(axis=1, skipna=False)

    # ------------------------------------------------------------------ 7
    banner("STEP 7: TARGET LABEL CREATION & DATA LEAKAGE PREVENTION")

    print(f"-> Defining 'target_toxo_pos' from SSTOXG >= {POS_CUTOFF} IU/mL...")
    merged['target_toxo_pos'] = (merged['SSTOXG'] >= POS_CUTOFF).astype(int)

    n_pos = int((merged['target_toxo_pos'] == 1).sum())
    n_neg = int((merged['target_toxo_pos'] == 0).sum())
    print(f"  Class 0 (Negative): {n_neg} ({n_neg / len(merged) * 100:.2f}%)")
    print(f"  Class 1 (Positive): {n_pos} ({n_pos / len(merged) * 100:.2f}%)")

    # Columns that must NOT enter the model:
    #   SSTOXG                              raw titer -> direct leakage of the label
    #   SMQ020, SMQ040                      replaced by smoking_status
    #   ALQ101, ALQ110, ALQ120Q, ALQ130     replaced by drinks_per_day
    #   DMDEDUC3                            folded into DMDEDUC2
    cols_to_drop = ['SSTOXG', 'SMQ020', 'SMQ040',
                    'ALQ101', 'ALQ110', 'ALQ120Q', 'ALQ130', 'DMDEDUC3']
    print(f"\n-> Dropping redundant / leakage columns: {cols_to_drop}")
    clean_df = merged.drop(columns=cols_to_drop)

    # ------------------------------------------------------------------ 8
    banner("STEP 8: FINAL VERIFICATION & EXPORT")

    print(f"Final dataset shape: {clean_df.shape} "
          f"(rows: {clean_df.shape[0]}, columns: {clean_df.shape[1]})")
    print(f"Columns ({len(clean_df.columns)}):")
    for i, col in enumerate(clean_df.columns, 1):
        print(f"  {i:2d}. {col:18} | dtype: {str(clean_df[col].dtype):8} | "
              f"missing: {clean_df[col].isna().sum():4d}")

    # Sanity check: the surviving values of each coded column should be only the
    # valid codes (no leftover 7/9/77/99).
    print("\nSurviving values in coded columns (should contain no Refused/DK codes):")
    for col in ['RIAGENDR', 'RIDRETH3', 'DMDBORN4', 'DMDEDUC2', 'HSD010',
                'SLQ050', 'smoking_status']:
        print(f"  {col:15}: {sorted(clean_df[col].dropna().unique().tolist())}")
    print(f"  {'SLD010H':15}: min={clean_df['SLD010H'].min()}, max={clean_df['SLD010H'].max()}")
    print(f"  {'drinks_per_day':15}: min={clean_df['drinks_per_day'].min()}, "
          f"max={clean_df['drinks_per_day'].max()}")
    print(f"  {'phq9_total_score':15}: min={clean_df['phq9_total_score'].min()}, "
          f"max={clean_df['phq9_total_score'].max()}")

    # Structural missingness by age band (helps spot questions skipped by age).
    print("\nShare of missing values per feature, by age band:")
    age_bins = pd.cut(clean_df['RIDAGEYR'], bins=[17, 19, 29, 39, 49, 59, 69, 120],
                      labels=['18-19', '20-29', '30-39', '40-49', '50-59', '60-69', '70+'])
    feature_cols = [c for c in clean_df.columns if c not in ('SEQN', 'target_toxo_pos')]
    miss_by_age = clean_df[feature_cols].isna().groupby(age_bins, observed=True).mean().round(3)
    with pd.option_context('display.max_columns', None, 'display.width', 200):
        print(miss_by_age)

    if WRITE_CSV:
        os.makedirs(PROCESSED_DIR, exist_ok=True)
        clean_df.to_csv(OUTPUT_CSV, index=False)
        print(f"\n[OK] Saved cleaned dataset to: {OUTPUT_CSV}")
    else:
        print("\n[NOTE] WRITE_CSV = False -> no file was written to disk.")
    print("  Reminder: drop 'SEQN' before training, and one-hot encode the nominal "
          "columns (RIAGENDR, RIDRETH3, DMDBORN4).")
    print("=" * 70)


if __name__ == '__main__':
    main()