"""Prepare analysis cohorts from the three database exports (read-only).

Steps: read exports -> score the 27-item FI-Lab (first values) -> select one admission per
patient (MIMIC-IV earliest; eICU random, R-compatible seed 20261022) -> outcome, covariates,
SOFA -> write private design matrices for the R main analysis and a cohort pickle for the
additional Python analyses.

Usage:
  python -I scripts/prepare.py --mimic mimic.xlsx --eicu eicu_max及min.xlsx \
      --eicu-first eicu首测值xlsx.xlsx --spec spec --output OUT
All outputs are patient-level and must stay private.
"""
import argparse, json, math, random
from pathlib import Path
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
for n in ['mimic', 'eicu', 'eicu-first', 'spec', 'output']: ap.add_argument('--' + n, required=True)
a = ap.parse_args()
OUT = Path(a.output); (OUT / 'data').mkdir(parents=True, exist_ok=True); (OUT / 'results').mkdir(parents=True, exist_ok=True)
spec = pd.read_csv(Path(a.spec) / 'item_crosswalk_and_thresholds.csv')

m = pd.read_excel(a.mimic); e = pd.read_excel(a.eicu); fe = pd.read_excel(a.eicu_first)
assert (e.patientunitstayid.values == fe.patientunitstayid.values).all(), 'eICU exports must share row order'

# ---------------- R-compatible random admission selection (R >= 3.6, Mersenne-Twister, Rejection) ----------------
class RRNG:
    """Reproduces set.seed(seed); sample.int(n, 1) from R."""
    def __init__(self, seed):
        s = seed & 0xffffffff
        for _ in range(50): s = (69069 * s + 1) & 0xffffffff
        st = []
        for _ in range(625): s = (69069 * s + 1) & 0xffffffff; st.append(s)
        self.g = random.Random(); self.g.setstate((3, tuple(st[1:]) + (624,), None))
    def unif(self):
        v = self.g.getrandbits(32) * 2.3283064365386963e-10
        if v <= 0: return 0.5 * 2.328306437080797e-10
        if 1 - v <= 0: return 1 - 0.5 * 2.328306437080797e-10
        return v
    def unif_index(self, dn):
        bits = int(math.ceil(math.log2(dn))) if dn > 1 else 0
        while True:
            v, n = 0, 0
            while n <= bits: v = 65536 * v + int(math.floor(self.unif() * 65536)); n += 16
            v &= (1 << bits) - 1
            if dn > v: return v

def select_eicu(e, seed=20261022):
    pid = e.uniquepid.astype(str).values; adm = e.patienthealthsystemstayid.values
    age = pd.to_numeric(e.age, errors='coerce').values
    elig = (~np.isnan(age)) & (age >= 18) & (e.unitvisitnumber.values == 1)
    cand = np.where(elig)[0]; ix = cand[np.lexsort((adm[cand], pid[cand]))]
    groups = {}
    for i in ix: groups.setdefault(pid[i], []).append(i)
    rng = RRNG(seed); keep = [g[rng.unif_index(len(g))] for _, g in sorted(groups.items())]
    sel = np.zeros(len(e), bool); sel[keep] = True
    return sel, elig

# ---------------- FI-Lab scoring ----------------
def flags(first, sex, co):
    F, V = {}, {}
    for r in spec.itertuples():
        b = r.mimic_column[:-4]
        if co == 'E': b = {'LDH': 'custom_乳酸脱氢酶', 'PaO2': 'pao2', 'PaCO2': 'paco2'}.get(r.item, b)
        x = pd.to_numeric(first[b + '_first'], errors='coerce').values
        if r.item == 'Hemoglobin': lo = np.where(sex == 'M', 14, 12.); hi = np.where(sex == 'M', 18, 16.)
        elif r.item == 'Creatinine': lo = np.where(sex == 'M', .74, .4); hi = np.where(sex == 'M', 1.35, 1.1)
        else: lo, hi = r.lower, r.upper
        F[r.item] = np.where(np.isnan(x), np.nan, ((x < lo) | (x > hi)).astype(float)); V[r.item] = x
    return pd.DataFrame(F, index=first.index), pd.DataFrame(V, index=first.index)

def eth(r):
    r = r.astype(str).str.upper()
    return np.select([r.str.contains('WHITE|CAUCASIAN'), r.str.contains('BLACK|AFRICAN'), r.str.contains('ASIAN')], ['White', 'Black', 'Asian'], 'Other/Unknown')

COV = ['hypertension', 'diabetes', 'congestive_heart_failure', 'severe_liver_disease', 'renal_disease', 'chronic_pulmonary_disease']
out, flow = {}, {}
# MIMIC-IV: earliest candidate admission per patient
d = m.copy(); fl = [('Candidate AMI ICU records', len(d))]
F, V = flags(d, d.gender.values, 'M')
d = d[(d.admission_age >= 18) & d.first_icu_stay.astype(str).str.lower().isin(['true', '1'])]; fl.append(('Adult, first ICU record per hospitalization', len(d)))
d = d.sort_values(['subject_id', 'admittime', 'icu_intime'], kind='mergesort').drop_duplicates('subject_id').sort_index(); fl.append(('One admission per patient', len(d)))
d['fol'] = (pd.to_datetime(d.dischtime) - pd.to_datetime(d.icu_intime)).dt.total_seconds() / 86400; st = d.hospital_expire_flag
d['age'] = d.admission_age; d['sexM'] = (d.gender == 'M') * 1.; d['eth'] = eth(d.race); d['id'] = d.stay_id; d['pid'] = d.subject_id
d['sev'] = pd.to_numeric(d.apsiii, errors='coerce'); d['cci_raw'] = d.charlson
out['M'] = (d, st, F, V, fl)
# eICU: one random candidate admission per patient
d = e.copy(); sex = d.gender.map({1: 'M', 0: 'F'}).values; fl = [('Candidate AMI ICU records', len(d))]
F, V = flags(fe, sex, 'E')
sel, elig = select_eicu(e); fl.append(('Adult, first ICU record per hospitalization', int(elig.sum())))
d = d[sel]; fl.append(('One admission per patient', len(d)))
d['fol'] = d.hospitaldischargeoffset / 1440; st = d.hosp_mort
d['cci_raw'] = d.charlson + 1.
se = pd.cut(d.age, [-np.inf, 49, 59, 69, 79, np.inf], labels=False); sm = pd.cut(d.age, [-np.inf, 50, 60, 70, 80, np.inf], labels=False)
d['charlson'] = d.charlson - se + sm + 1; d['sexM'] = (d.gender == 1) * 1.; d['eth'] = eth(d.ethnicity); d['id'] = d.patientunitstayid; d['pid'] = d.uniquepid
d['sev'] = pd.to_numeric(d.apacheIV, errors='coerce')
out['E'] = (d, st, F, V, fl)

res, designs = {}, {}
for co, (d, st, F, V, fl) in out.items():
    d = d.copy()
    d = d[st.loc[d.index].isin([0, 1]) & d.fol.notna() & (d.fol > 0)]; fl.append(('Valid hospital outcome and timing', len(d)))
    d['death30'] = ((st.loc[d.index] == 1) & (d.fol <= 30)) * 1.
    d = d[d[COV + ['charlson', 'age']].notna().all(1)]; fl.append(('Complete covariates', len(d)))
    Fi = F.loc[d.index]; d['ncomp'] = Fi.notna().sum(1).values; d['nabn'] = Fi.sum(1).values
    d['fi_any'] = np.where(d.ncomp > 0, d.nabn / np.maximum(d.ncomp, 1), np.nan)
    d['sofa'] = pd.to_numeric(d.sofa_dynamic_max, errors='coerce')
    h = pd.to_numeric(d.height, errors='coerce'); w = pd.to_numeric(d.weight, errors='coerce'); bmi = w / (h / 100) ** 2
    bmi[(h < 100) | (h > 250) | (w < 20) | (w > 400) | (bmi < 10) | (bmi > 80)] = np.nan
    d['bmi30'] = (bmi > 30) * 1.; d['bmiM'] = bmi.isna() * 1.
    for k in ['Black', 'Asian', 'Other/Unknown']: d['eth_' + k] = (d.eth == k) * 1.
    cr = pd.to_numeric(d['creatinine_max'], errors='coerce')
    crs = pd.cut(cr, [-np.inf, 1.2, 2.0, 3.5, 5.0, np.inf], right=False, labels=[0, 1, 2, 3, 4]).astype(float).fillna(0)
    d['sofa_cr'] = d.sofa - pd.to_numeric(d.sofa_dynamic_renal_max, errors='coerce') + crs
    fl.append(('FI-Lab >=18 components', int((d.ncomp >= 18).sum())))
    fl.append(('FI-Lab >=18 components and SOFA available', int(((d.ncomp >= 18) & d.sofa.notna()).sum())))
    k = d[(d.ncomp >= 18) & d.sofa.notna()]
    fl.append(('Excluded within 24 hours: deaths', int(((k.fol <= 1) & (k.death30 == 1)).sum())))
    fl.append(('Excluded within 24 hours: other discharges', int(((k.fol <= 1) & (k.death30 == 0)).sum())))
    an = k[k.fol > 1]; fl.append(('Analysis cohort', len(an))); fl.append(('Analysis cohort deaths', int(an.death30.sum())))
    res[co] = dict(d=d, F=Fi, V=V.loc[d.index], flow=fl)
    X = pd.DataFrame({'const': 1.0, 'age': k.age, 'sexM': k.sexM, 'ethBlack': k['eth_Black'], 'ethAsian': k['eth_Asian'], 'ethOther/Unknown': k['eth_Other/Unknown'],
                      'bmi>30': k.bmi30, 'bmiMissing': k.bmiM, **{c: k[c] * 1. for c in COV}, 'charlson': k.charlson, 'sofa': k.sofa,
                      'fi01': k.fi_any * 10, 'death30': k.death30, 'followup_days': k.fol, 'id': k.id, 'pid': k.pid})
    cols = pd.read_csv(Path(a.spec) / 'design_columns.csv', nrows=0).columns.tolist(); assert list(X.columns) == cols
    X.to_csv(OUT / 'data' / f'{co}_design.csv', index=False)
pd.to_pickle(res, OUT / 'data' / 'cohorts.pkl')
json.dump({co: res[co]['flow'] for co in res}, open(OUT / 'results' / 'cohort_flow.json', 'w'), indent=1)
for co in res: print(co, res[co]['flow'])
