"""Export the design columns needed by threshold.R (all patients, before the >=18-component rule).
Usage: python -I threshold_export.py COHORTS_PKL OUTDIR   -> OUTDIR/{M,E}_all.csv (patient-level; keep private)"""
import sys
import pandas as pd
pk, od = sys.argv[1:3]
r = pd.read_pickle(pk)
cols = ['age', 'sexM', 'eth_Black', 'eth_Asian', 'eth_Other/Unknown', 'bmi30', 'bmiM', 'hypertension', 'diabetes', 'congestive_heart_failure',
        'severe_liver_disease', 'renal_disease', 'chronic_pulmonary_disease', 'charlson', 'sofa', 'sofa_cr', 'fi_any', 'ncomp', 'death30', 'fol']
for co in ['M', 'E']:
    r[co]['d'][cols].to_csv(f'{od}/{co}_all.csv', index=False)
print('ok')
