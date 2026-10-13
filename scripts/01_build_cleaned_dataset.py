"""
EXACT WORKFLOW:
1. Load 7 official CDC NHANES .xpt files from data/raw/NHANES/
2. Find the exact participant intersection on 'SEQN' 
3. Extract only clinically justified features (Demographics, Health, Sleep, Mood, Lifestyle)
4. Clean special CDC codes (7/77/777=Refused, 9/99/999=Don't Know) -> NaN
5. Clean SAS float artifacts (5.3976e-79) -> 0.0
6. Resolve survey skip logic:
   - Smoking (SMQ020 + SMQ040 -> smoking_status: 0=Never, 1=Former, 2=Current)
   - Alcohol (ALQ101 + ALQ130 -> drinks_per_day: non-drinkers = 0.0)
   - Education for 18-19 year olds (DMDEDUC2 -> High School equivalent)
7. Compute composite clinical score:
   - phq9_total_score (Sum of 9 depression items, scale 0 to 27)
8. Define binary target variable:
   - target_toxo_pos = 1 if SSTOXG >= 33.0 IU/mL, else 0
9. Drop raw lab titer (SSTOXG) to prevent 100% data leakage
10. Save final dataset to: data/processed/nhanes_adult_toxo_clean.csv
================================================================================
"""

import os
import numpy as np
import pandas as pd

def main():
    raw_dir = os.path.join('data', 'raw', 'NHANES')
    processed_dir = os.path.join('data', 'processed')
    os.makedirs(processed_dir, exist_ok=True)
    output_csv = os.path.join(processed_dir, 'nhanes_adult_toxo_clean.csv')

    print("=" * 70)
    print("STEP 1: LOADING RAW CDC NHANES FILES & EXTRACTING SPECIFIC COLUMNS")
    print("=" * 70)

    # 1. Target: Toxoplasma IgG antibodies
    print("-> Reading SSTOXO_H.xpt (Toxoplasma IgG)...")
    df_toxo = pd.read_sas(os.path.join(raw_dir, 'SSTOXO_H.xpt'), format='xport')[['SEQN', 'SSTOXG']]

    # 2. Demographics
    print("-> Reading DEMO_H.xpt (Demographics)...")
    demo_cols = ['SEQN', 'RIDAGEYR', 'RIAGENDR', 'RIDRETH3', 'DMDBORN4', 'DMDEDUC2', 'INDFMPIR']
    df_demo = pd.read_sas(os.path.join(raw_dir, 'DEMO_H.xpt'), format='xport')[demo_cols]

    # 3. General Health
    print("-> Reading HSQ_H.xpt (General Health)...")
    df_hsq = pd.read_sas(os.path.join(raw_dir, 'HSQ_H.xpt'), format='xport')[['SEQN', 'HSD010']]

    # 4. Sleep Disorders
    print("-> Reading SLQ_H.xpt (Sleep Disorders)...")
    df_slq = pd.read_sas(os.path.join(raw_dir, 'SLQ_H.xpt'), format='xport')[['SEQN', 'SLD010H', 'SLQ050']]

    # 5. Depression PHQ-9 Screener
    print("-> Reading DPQ_H.xpt (Depression Screener)...")
    dpq_cols = ['SEQN', 'DPQ010', 'DPQ020', 'DPQ030', 'DPQ040', 'DPQ050', 'DPQ060', 'DPQ070', 'DPQ080', 'DPQ090']
    df_dpq = pd.read_sas(os.path.join(raw_dir, 'DPQ_H.xpt'), format='xport')[dpq_cols]

    # 6. Smoking Questionnaire
    print("-> Reading SMQ_H.xpt (Smoking Questionnaire)...")
    df_smq = pd.read_sas(os.path.join(raw_dir, 'SMQ_H.xpt'), format='xport')[['SEQN', 'SMQ020', 'SMQ040']]

    # 7. Alcohol Questionnaire
    print("-> Reading ALQ_H.xpt (Alcohol Questionnaire)...")
    df_alq = pd.read_sas(os.path.join(raw_dir, 'ALQ_H.xpt'), format='xport')[['SEQN', 'ALQ101', 'ALQ130']]

    print("\n" + "=" * 70)
    print("STEP 2: FINDING THE ADULT COMMON INTERSECTION ON 'SEQN'")
    print("=" * 70)

    # Calculate exact intersection of SEQN across all 7 files
    seqn_sets = [
        set(df_toxo['SEQN']),
        set(df_demo['SEQN']),
        set(df_hsq['SEQN']),
        set(df_slq['SEQN']),
        set(df_dpq['SEQN']),
        set(df_smq['SEQN']),
        set(df_alq['SEQN'])
    ]
    common_seqn = set.intersection(*seqn_sets)
    print(f"Total participants with complete overlap across all 7 files: {len(common_seqn)}")

    #join
    merged = df_toxo[df_toxo['SEQN'].isin(common_seqn)].copy()
    for df_next in [df_demo, df_hsq, df_slq, df_dpq, df_smq, df_alq]:
        merged = merged.merge(df_next, on='SEQN', how='inner')

    print(f"Merged cohort shape: {merged.shape} (Rows: {merged.shape[0]}, Columns: {merged.shape[1]})")

    print("\n" + "=" * 70)
    print("STEP 3: CLEANING SAS FLOAT ZEROES (5.3976e-79 -> 0.0)")
    print("=" * 70)

    sas_zero = 5.397605346934028e-79
    for col in dpq_cols[1:] + ['INDFMPIR']:
        mask = np.isclose(merged[col].fillna(-999), sas_zero)
        count = mask.sum()
        if count > 0:
            merged.loc[mask, col] = 0.0
            print(f"  Fixed {count:4d} SAS zeroes in column: {col}")

    print("\n" + "=" * 70)
    print("STEP 4: CLEANING CDC SPECIAL CODES (7=Refused, 9=Don't Know -> NaN)")
    print("=" * 70)

    # For categorical/ordinal flags, 7/77/777 = Refused, 9/99/999 = Don't Know
    # We convert these to NaN so ML algorithms don't treat them as real high numbers!
    categorical_cols_to_clean = {
        'DMDBORN4': [7, 9],
        'DMDEDUC2': [7, 9],
        'HSD010':   [7, 9],
        'SLQ050':   [7, 9],
        'SMQ020':   [7, 9],
        'ALQ101':   [7, 9],
    }
    # Add all 9 PHQ-9 items
    for col in dpq_cols[1:]:
        categorical_cols_to_clean[col] = [7, 9]

    # For Sleep Hours (SLD010H), 99 = Don't know (7 is 7 hours of sleep, so DO NOT remove 7!)
    categorical_cols_to_clean['SLD010H'] = [99]

    for col, special_codes in categorical_cols_to_clean.items():
        mask = merged[col].isin(special_codes)
        count = mask.sum()
        if count > 0:
            merged.loc[mask, col] = np.nan
            print(f"  Converted {count:3d} Refused/Don't-know entries to NaN in: {col}")

    print("\n" + "=" * 70)
    print("STEP 5: RESOLVING SURVEY SKIP LOGIC & FEATURE ENGINEERING")
    print("=" * 70)

    # A. SMOKING STATUS
    # SMQ020: 1=Yes (smoked 100 cigs), 2=No
    # SMQ040: 1=Every day, 2=Some days, 3=Not at all (only asked if SMQ020==1)
    print("-> Engineering 'smoking_status'...")
    # 0 = Never smoker, 1 = Former smoker, 2 = Current smoker
    smoking_status = pd.Series(np.nan, index=merged.index)
    smoking_status[merged['SMQ020'] == 2] = 0.0  # Never smoker
    smoking_status[(merged['SMQ020'] == 1) & (merged['SMQ040'] == 3)] = 1.0  # Former smoker
    smoking_status[(merged['SMQ020'] == 1) & (merged['SMQ040'].isin([1, 2]))] = 2.0  # Current smoker
    merged['smoking_status'] = smoking_status

    # B. ALCOHOL DRINKS PER DAY
    # ALQ101: 1=Yes (had 12 drinks/year), 2=No
    # ALQ130: Average drinks per day on days drank
    print("-> Engineering 'drinks_per_day'...")
    drinks_per_day = merged['ALQ130'].copy()
    drinks_per_day[merged['ALQ101'] == 2] = 0.0  # Non-drinkers drink 0 drinks per day
    # Clean refusal codes (777, 999) if any
    drinks_per_day[drinks_per_day.isin([77, 99, 777, 999])] = np.nan
    merged['drinks_per_day'] = drinks_per_day

    # C. EDUCATION LEVEL FOR 18-19 YEAR OLDS
    # DMDEDUC2 is NaN for all 252 18-19 year olds because it's only asked for 20+
    print("-> Resolving education for 18-19 year olds (DMDEDUC2)...")
    # Code 3 in DMDEDUC2 = High school graduate / GED / High school student
    merged.loc[merged['DMDEDUC2'].isna() & merged['RIDAGEYR'].isin([18, 19]), 'DMDEDUC2'] = 3.0

    # D. COMPOSITE PHQ-9 DEPRESSION TOTAL SCORE
    print("-> Engineering 'phq9_total_score' (Sum of 9 items, 0 to 27)...")
    dpq_items = dpq_cols[1:]
    merged['phq9_total_score'] = merged[dpq_items].sum(axis=1, skipna=False)

    print("\n" + "=" * 70)
    print("STEP 6: TARGET LABEL CREATION & DATA LEAKAGE PREVENTION")
    print("=" * 70)

    # CDC Cutoff: SSTOXG >= 33.0 IU/mL is Positive (1), else Negative (0)
    print("-> Defining 'target_toxo_pos' from SSTOXG >= 33.0 IU/mL...")
    merged['target_toxo_pos'] = (merged['SSTOXG'] >= 33.0).astype(int)

    n_pos = int((merged['target_toxo_pos'] == 1).sum())
    n_neg = int((merged['target_toxo_pos'] == 0).sum())
    print(f"  Class 0 (Negative): {n_neg} ({n_neg/len(merged)*100:.2f}%)")
    print(f"  Class 1 (Positive): {n_pos} ({n_pos/len(merged)*100:.2f}%)")

    # DROPPING COLUMNS THAT MUST NOT ENTER THE MACHINE LEARNING MODEL:
    # 1. SSTOXG: Raw antibody titer must be dropped to prevent 100% data leakage!
    # 2. SMQ020, SMQ040: Raw smoking columns replaced by 'smoking_status'
    # 3. ALQ101, ALQ130: Raw alcohol columns replaced by 'drinks_per_day'
    cols_to_drop = ['SSTOXG', 'SMQ020', 'SMQ040', 'ALQ101', 'ALQ130']
    print(f"\n-> Dropping redundant / leakage columns: {cols_to_drop}")
    clean_df = merged.drop(columns=cols_to_drop)

    print("\n" + "=" * 70)
    print("STEP 7: FINAL VERIFICATION & EXPORT TO CSV")
    print("=" * 70)

    print(f"Final Clean Dataset Shape: {clean_df.shape} (Rows: {clean_df.shape[0]}, Features: {clean_df.shape[1]})")
    print(f"Columns in clean dataset ({len(clean_df.columns)}):")
    for i, col in enumerate(clean_df.columns, 1):
        print(f"  {i:2d}. {col:18} | Type: {str(clean_df[col].dtype):8} | Missing (NaN): {clean_df[col].isna().sum():4d}")

    # Export to CSV (PAUSED per user request until results are verified)
    # clean_df.to_csv(output_csv, index=False)
    print("\n[NOTE] CSV generation is paused per your request. No file was written to disk.")
    print("=" * 70)

if __name__ == '__main__':
    main()

