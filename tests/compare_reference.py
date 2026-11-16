"""Compare a run against reference_results (aggregate results only).
Usage: python -I tests/compare_reference.py OUTPUT_DIR"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[1]; OUT = Path(sys.argv[1]); ok = True
def cmp_csv(ref, new, keys, tol=1e-6):
    global ok
    if not new.exists(): print('SKIP (not produced):', new.name); return
    a, b = pd.read_csv(ref), pd.read_csv(new)
    num = [c for c in a.columns if c not in keys and pd.api.types.is_numeric_dtype(a[c])]
    m = a.merge(b, on=keys, suffixes=('_ref', '_new'))
    assert len(m) == len(a), f'{new.name}: row mismatch'
    diff = max((m[c + '_ref'] - m[c + '_new']).abs().max() for c in num)
    flag = diff <= tol; ok &= flag; print(('OK  ' if flag else 'DIFF'), new.name, f'max |diff| = {diff:.2e}')
def cmp_json(ref, new, tol=1e-6):
    global ok
    def flat(x, k=''):
        if isinstance(x, dict):
            for kk, v in x.items(): yield from flat(v, k + '/' + kk)
        elif isinstance(x, list):
            for i, v in enumerate(x): yield from flat(v, f'{k}[{i}]')
        elif isinstance(x, (int, float)): yield k, float(x)
    a, b = dict(flat(json.load(open(ref)))), dict(flat(json.load(open(new))))
    diff = max(abs(a[k] - b[k]) for k in a); flag = diff <= tol; ok &= flag
    print(('OK  ' if flag else 'DIFF'), new.name, f'max |diff| = {diff:.2e}')
R, N = ROOT / 'reference_results', OUT / 'results'
cmp_csv(R / 'main/external_validation.csv', N / 'main/external_validation.csv', ['scenario', 'metric'])
cmp_csv(R / 'main/internal_validation.csv', N / 'main/internal_validation.csv', ['scenario', 'metric'])
cmp_csv(R / 'main/exclusion_audit.csv', N / 'main/exclusion_audit.csv', ['cohort'])
cmp_csv(R / 'additional/external_sensitivity.csv', N / 'additional/external_sensitivity.csv', ['analysis'])
cmp_csv(R / 'additional/renal_sofa.csv', N / 'additional/renal_sofa.csv', ['analysis'])
cmp_csv(R / 'additional/components.csv', N / 'additional/components.csv', ['cohort', 'item'])
cmp_csv(R / 'additional/threshold_performance.csv', N / 'additional/threshold_performance.csv', ['threshold', 'sofa', 'metric'])
cmp_csv(R / 'additional/included_excluded.csv', N / 'additional/included_excluded.csv', ['cohort', 'group', 'variable'])
cmp_csv(R / 'additional/risk_grid.csv', N / 'additional/risk_grid.csv', ['cohort', 'sofa', 'fi_tertile'])
cmp_csv(R / 'additional/within_sofa_or.csv', N / 'additional/within_sofa_or.csv', ['cohort', 'sofa'])
cmp_csv(R / 'additional/eicu_region_validation.csv', N / 'additional/eicu_region_validation.csv', ['region'])
for f in ['within_severity.json', 'associations.json', 'overlap.json', 'clinical.json', 'heterogeneity.json']: cmp_json(R / 'additional' / f, N / 'additional' / f)
fl_ref, fl_new = json.load(open(R / 'cohort_flow.json')), json.load(open(N / 'cohort_flow.json'))
flag = fl_ref == fl_new; ok &= flag; print(('OK  ' if flag else 'DIFF'), 'cohort_flow.json')
sys.exit(0 if ok else 1)
