"""Clinical translation analyses (NumPy only):
(1) observed mortality by SOFA category x FI-Lab tertile (tertile cut-points from MIMIC-IV development);
(2) classification at a 10% predicted-risk threshold, frozen models in eICU, with paired bootstrap;
(3) standardized predicted risk at SOFA 6 with FI-Lab 0.3 vs 0.6 (frozen FI-Lab model, eICU covariates);
(4) three-variable model (age, SOFA, FI-Lab) developed in MIMIC-IV and applied unchanged to eICU.
Usage: python -I clinical.py OUT/data/cohorts.pkl OUT/results/additional"""
import sys, json, math, numpy as np, pandas as pd
PK, OD = sys.argv[1:3]; res = pd.read_pickle(PK); B = 1000

def cohort(co):
    d = res[co]['d']; return d[(d.ncomp >= 18) & d.sofa.notna() & (d.fol > 1)].copy()
COV = ['hypertension', 'diabetes', 'congestive_heart_failure', 'severe_liver_disease', 'renal_disease', 'chronic_pulmonary_disease']
def design(d, fi=True, sofa=None, fiv=None):
    s = d.sofa.values if sofa is None else np.full(len(d), sofa, float)
    cols = [np.ones(len(d)), d.age, d.sexM, d['eth_Black'], d['eth_Asian'], d['eth_Other/Unknown'], d.bmi30, d.bmiM] + [d[c] * 1. for c in COV] + [d.charlson, s]
    if fi: cols.append(d.fi_any.values * 10 if fiv is None else np.full(len(d), fiv * 10))
    return np.column_stack(cols).astype(float)
def logit(X, y):
    b = np.zeros(X.shape[1])
    for _ in range(100):
        p = 1 / (1 + np.exp(-X @ b)); s = np.linalg.solve((X * (p * (1 - p))[:, None]).T @ X, X.T @ (y - p)); b += s
        if np.abs(s).max() < 1e-10: break
    return b
def auc(y, p):
    r = pd.Series(p).rank().values; n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n0 * n1)
def wilson(k, n):
    z = 1.959964; p = k / n; c = (p + z * z / (2 * n)) / (1 + z * z / n); h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h
P = lambda X, b: 1 / (1 + np.exp(-X @ b))

M, E = cohort('M'), cohort('E'); ym, ye = M.death30.values, E.death30.values
b0 = logit(design(M, fi=False), ym); b1 = logit(design(M), ym)
p0, p1 = P(design(E, fi=False), b0), P(design(E), b1)
out = {}

# (1) risk grid
cuts = list(np.quantile(M.fi_any, [1 / 3, 2 / 3]))
def sofa_cat(s): return np.select([s <= 3, s <= 7], ['SOFA 0–3', 'SOFA 4–7'], 'SOFA ≥8')
def fi_cat(f): return np.select([f <= cuts[0], f <= cuts[1]], ['Low', 'Middle'], 'High')
grid = []
for co, d in [('MIMIC-IV', M), ('eICU', E)]:
    sc, fc = sofa_cat(d.sofa.values), fi_cat(d.fi_any.values)
    for s in ['SOFA 0–3', 'SOFA 4–7', 'SOFA ≥8']:
        for f in ['Low', 'Middle', 'High']:
            m = (sc == s) & (fc == f); n = int(m.sum()); k = int(d.death30.values[m].sum())
            lo, hi = wilson(k, n) if n else (np.nan, np.nan)
            grid.append(dict(cohort=co, sofa=s, fi_tertile=f, n=n, deaths=k, mortality=k / n if n else np.nan, lo=lo, hi=hi))
pd.DataFrame(grid).to_csv(OD + '/risk_grid.csv', index=False)
out['fi_tertile_cuts'] = cuts

# (2) 10% threshold
def cls(y, p, t=0.10):
    pos = p >= t; tp = (pos & (y == 1)).sum(); fp = (pos & (y == 0)).sum(); fn = (~pos & (y == 1)).sum(); tn = (~pos & (y == 0)).sum()
    return dict(sens=tp / (tp + fn), spec=tn / (tn + fp), ppv=tp / (tp + fp), npv=tn / (tn + fn), flagged=pos.mean())
c0, c1 = cls(ye, p0), cls(ye, p1); rate = ye.mean()
def per1000(a, b): return dict(extra_deaths_identified=(b['sens'] - a['sens']) * rate * 1000, change_false_positives=-(b['spec'] - a['spec']) * (1 - rate) * 1000)
rng = np.random.default_rng(20261030); bs = []
for _ in range(B):
    i = rng.integers(0, len(ye), len(ye)); a, b = cls(ye[i], p0[i]), cls(ye[i], p1[i]); r_ = ye[i].mean()
    bs.append([b['sens'] - a['sens'], b['spec'] - a['spec'], (b['sens'] - a['sens']) * r_ * 1000, -(b['spec'] - a['spec']) * (1 - r_) * 1000])
bs = np.array(bs); q = np.percentile(bs, [2.5, 97.5], axis=0)
out['threshold10'] = dict(without=c0, with_=c1, event_rate=rate, **per1000(c0, c1),
                          d_sens_ci=q[:, 0].tolist(), d_spec_ci=q[:, 1].tolist(), extra_deaths_ci=q[:, 2].tolist(), change_fp_ci=q[:, 3].tolist(),
                          d_sens=c1['sens'] - c0['sens'], d_spec=c1['spec'] - c0['spec'])

# (3) standardized risk at SOFA 6, FI-Lab 0.3 vs 0.6
r3 = P(design(E, sofa=6, fiv=0.3), b1).mean(); r6 = P(design(E, sofa=6, fiv=0.6), b1).mean()
out['standardized_sofa6'] = dict(fi_0_3=r3, fi_0_6=r6, difference=r6 - r3)

# (4) three-variable model
def X3(d, fi=True): return np.column_stack([np.ones(len(d)), d.age, d.sofa] + ([d.fi_any * 10] if fi else [])).astype(float)
a0, a1 = logit(X3(M, False), ym), logit(X3(M), ym)
q0, q1 = P(X3(E, False), a0), P(X3(E), a1)
def met(y, p):
    lp = np.log(p / (1 - p)); s = logit(np.column_stack([np.ones(len(y)), lp]), y)[1]
    return dict(AUC=auc(y, p), OE=y.sum() / p.sum(), slope=s, Brier=np.mean((p - y) ** 2))
m30, m31, mfull = met(ye, q0), met(ye, q1), met(ye, p1)
rng = np.random.default_rng(20261031); bs = []
for _ in range(B):
    i = rng.integers(0, len(ye), len(ye)); y = ye[i]
    bs.append([auc(y, q0[i]), auc(y, q1[i]), auc(y, q1[i]) - auc(y, q0[i]), auc(y, p1[i]) - auc(y, q1[i])])
bs = np.array(bs); q = np.percentile(bs, [2.5, 97.5], axis=0)
out['three_variable'] = dict(coef_with=dict(zip(['intercept', 'age', 'sofa', 'fi01'], a1.tolist())), age_sofa=m30, age_sofa_fi=m31, full=mfull,
                             auc_age_sofa_ci=q[:, 0].tolist(), auc_age_sofa_fi_ci=q[:, 1].tolist(), dauc=m31['AUC'] - m30['AUC'], dauc_ci=q[:, 2].tolist(),
                             full_minus_3var=mfull['AUC'] - m31['AUC'], full_minus_3var_ci=q[:, 3].tolist())
json.dump(out, open(OD + '/clinical.json', 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float)); print(pd.DataFrame(grid).round(3).to_string())
