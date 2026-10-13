"""
================================================================================
SCRIPT: check_toxo_positives_negatives.py
PURPOSE: Check and verify the exact number of Toxoplasma Positives and Negatives
         in both the full raw NHANES file (N = 6,578) and the adult cohort (N = 5,006).
================================================================================
"""

import os
import pandas as pd

def main():
    raw_dir = os.path.join('data', 'raw', 'NHANES')
    toxo_file = os.path.join(raw_dir, 'SSTOXO_H.xpt')
    demo_file = os.path.join(raw_dir, 'DEMO_H.xpt')

    print("=" * 70)
    print("1. FULL RAW TOXOPLASMA DATASET (SSTOXO_H.xpt)")
    print("=" * 70)

    df_toxo_full = pd.read_sas(toxo_file, format='xport')
    total_tested = len(df_toxo_full)
    
    # Official CDC Cutoffs
    raw_pos = (df_toxo_full['SSTOXG'] >= 33.0).sum()
    raw_equiv = ((df_toxo_full['SSTOXG'] >= 27.0) & (df_toxo_full['SSTOXG'] < 33.0)).sum()
    raw_neg = (df_toxo_full['SSTOXG'] < 27.0).sum()
    raw_missing = df_toxo_full['SSTOXG'].isna().sum()

    print(f"Total participants tested (all ages 6+): {total_tested:,}")
    print(f"  - Positive  (IgG >= 33.0 IU/mL):      {raw_pos:5,d} ({raw_pos/total_tested*100:5.2f}%)")
    print(f"  - Equivocal (27.0 <= IgG < 33.0):     {raw_equiv:5,d} ({raw_equiv/total_tested*100:5.2f}%)")
    print(f"  - Negative  (IgG < 27.0 IU/mL):       {raw_neg:5,d} ({raw_neg/total_tested*100:5.2f}%)")
    print(f"  - Missing   (No test result):         {raw_missing:5,d}")

    print("\n" + "=" * 70)
    print("2. ADULT COHORT OVERLAP (N = 5,006)")
    print("=" * 70)

    # Load Demographics to get age, gender, country of birth
    df_demo = pd.read_sas(demo_file, format='xport')[['SEQN', 'RIDAGEYR', 'RIAGENDR', 'DMDBORN4']]

    # Find common SEQN across all 7 questionnaire files
    files = ['SSTOXO_H.xpt', 'DEMO_H.xpt', 'HSQ_H.xpt', 'SLQ_H.xpt', 'DPQ_H.xpt', 'SMQ_H.xpt', 'ALQ_H.xpt']
    seqn_sets = [set(pd.read_sas(os.path.join(raw_dir, f), format='xport')['SEQN']) for f in files]
    common_seqn = set.intersection(*seqn_sets)

    # Merge cohort
    cohort = df_toxo_full[df_toxo_full['SEQN'].isin(common_seqn)].merge(df_demo, on='SEQN', how='inner').copy()
    cohort['target_toxo_pos'] = (cohort['SSTOXG'] >= 33.0).astype(int)

    n_cohort = len(cohort)
    pos_cohort = int((cohort['target_toxo_pos'] == 1).sum())
    neg_cohort = int((cohort['target_toxo_pos'] == 0).sum())

    print(f"Total adults in our final study cohort:   {n_cohort:,}")
    print(f"  - Positives (Class 1 - Infected):       {pos_cohort:5,d} ({pos_cohort/n_cohort*100:5.2f}%)")
    print(f"  - Negatives (Class 0 - Uninfected):     {neg_cohort:5,d} ({neg_cohort/n_cohort*100:5.2f}%)")

    print("\n" + "=" * 70)
    print("3. DEMOGRAPHIC BREAKDOWN OF POSITIVE CASES IN OUR ADULT COHORT")
    print("=" * 70)

    # Gender breakdown
    pos_male = int(((cohort['target_toxo_pos'] == 1) & (cohort['RIAGENDR'] == 1)).sum())
    pos_female = int(((cohort['target_toxo_pos'] == 1) & (cohort['RIAGENDR'] == 2)).sum())
    total_male = int((cohort['RIAGENDR'] == 1).sum())
    total_female = int((cohort['RIAGENDR'] == 2).sum())

    print("A. By Gender:")
    print(f"   - Males:   {pos_male} positives out of {total_male} ({pos_male/total_male*100:.2f}% positive rate)")
    print(f"   - Females: {pos_female} positives out of {total_female} ({pos_female/total_female*100:.2f}% positive rate)")

    # Age breakdown
    cohort['age_group'] = pd.cut(
        cohort['RIDAGEYR'],
        bins=[17, 35, 55, 100],
        labels=['Young Adults (18-35)', 'Middle-aged (36-55)', 'Older Adults (56+)']
    )
    print("\nB. By Age Group:")
    for grp in ['Young Adults (18-35)', 'Middle-aged (36-55)', 'Older Adults (56+)']:
        sub = cohort[cohort['age_group'] == grp]
        p = int((sub['target_toxo_pos'] == 1).sum())
        t = len(sub)
        print(f"   - {grp:22}: {p:3d} positives out of {t:5d} ({p/t*100:5.2f}% positive rate)")

    # Country of birth breakdown
    print("\nC. By Country of Birth (DMDBORN4):")
    # 1 = Born in 50 US States, 2 = Others (Foreign Born)
    us_sub = cohort[cohort['DMDBORN4'] == 1]
    foreign_sub = cohort[cohort['DMDBORN4'] == 2]
    p_us = int((us_sub['target_toxo_pos'] == 1).sum())
    p_foreign = int((foreign_sub['target_toxo_pos'] == 1).sum())
    print(f"   - US Born:       {p_us:3d} positives out of {len(us_sub):5d} ({p_us/len(us_sub)*100:5.2f}% positive rate)")
    print(f"   - Foreign Born:  {p_foreign:3d} positives out of {len(foreign_sub):5d} ({p_foreign/len(foreign_sub)*100:5.2f}% positive rate)")

    print("=" * 70)

if __name__ == '__main__':
    main()

