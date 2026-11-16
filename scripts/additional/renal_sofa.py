"""Table 2 (lower block): frozen-model performance with the SOFA renal subscore recalculated from
maximum creatinine (-6 to <24 h) in both databases. NumPy only.
Usage: python -I renal_sofa.py OUT/data/cohorts.pkl OUT/results/additional"""
import sys, math, numpy as np, pandas as pd
PK, OD = sys.argv[1:3]; res = pd.read_pickle(PK)

def cohort(co):
    d = res[co]['d']; return d[(d.ncomp >= 18) & d.sofa.notna() & (d.fol > 1)]
def design(d, sofa):
    return np.column_stack([np.ones(len(d)), d.age, d.sexM, d['eth_Black'], d['eth_Asian'], d['eth_Other/Unknown'], d.bmi30, d.bmiM]
                           + [d[c] * 1. for c in ['hypertension', 'diabetes', 'congestive_heart_failure', 'severe_liver_disease', 'renal_disease', 'chronic_pulmonary_disease']]
                           + [d.charlson, d[sofa]]).astype(float)
def logit(X, y):
    b = np.zeros(X.shape[1])
    for _ in range(60):
        p = 1 / (1 + np.exp(-X @ b)); g = X.T @ (y - p); H = (X * (p * (1 - p))[:, None]).T @ X
        s = np.linalg.solve(H, g); b += s
        if np.abs(s).max() < 1e-12: break
    return b
def auc(y, p):
    r = pd.Series(p).rank().values; n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n0 * n1)
def met(y, lp):
    p = 1 / (1 + np.exp(-lp)); a = 0.
    for _ in range(60):
        q = 1 / (1 + np.exp(-(lp + a))); a += (y - q).sum() / (q * (1 - q)).sum()
    s = logit(np.column_stack([np.ones(len(y)), lp]), y)[1]
    return dict(AUC=auc(y, p), Brier=np.mean((p - y) ** 2), OE=y.sum() / p.sum(), Intercept=a, Slope=s)
def delong(y, p1, p2):
    pos, neg = y == 1, y == 0
    def comp(p):
        x, z = p[pos], p[neg]
        V10 = np.array([((xi > z).sum() + 0.5 * (xi == z).sum()) / len(z) for xi in x])
        V01 = np.array([((x > zj).sum() + 0.5 * (x == zj).sum()) / len(x) for zj in z])
        return V10, V01
    v10a, v01a = comp(p1); v10b, v01b = comp(p2)
    S10 = np.cov(np.vstack([v10a, v10b])); S01 = np.cov(np.vstack([v01a, v01b]))
    var = (S10[0, 0] + S10[1, 1] - 2 * S10[0, 1]) / len(v10a) + (S01[0, 0] + S01[1, 1] - 2 * S01[0, 1]) / len(v01a)
    z = (v10a.mean() - v10b.mean()) / math.sqrt(var); return math.erfc(abs(z) / math.sqrt(2))

M, E = cohort('M'), cohort('E'); ym, ye = M.death30.values, E.death30.values
out = []
for label, sofa in [('Exported SOFA (main)', 'sofa'), ('Creatinine-based renal SOFA', 'sofa_cr')]:
    Xm, Xe = design(M, sofa), design(E, sofa); fm, fe = M.fi_any.values * 10, E.fi_any.values * 10
    b0 = logit(Xm, ym); b1 = logit(np.column_stack([Xm, fm]), ym)
    lp0 = Xe @ b0; lp1 = np.column_stack([Xe, fe]) @ b1
    m0, m1 = met(ye, lp0), met(ye, lp1)
    rng = np.random.default_rng(20261029); B = {k: [] for k in ['AUC0', 'AUC1', 'dAUC', 'Br0', 'Br1', 'dBr', 'OE0', 'OE1', 'I0', 'I1', 'S0', 'S1']}
    for _ in range(1000):
        i = rng.integers(0, len(ye), len(ye)); a, b = met(ye[i], lp0[i]), met(ye[i], lp1[i])
        for k, v in zip(B, [a['AUC'], b['AUC'], b['AUC'] - a['AUC'], a['Brier'], b['Brier'], a['Brier'] - b['Brier'], a['OE'], b['OE'], a['Intercept'], b['Intercept'], a['Slope'], b['Slope']]): B[k].append(v)
    row = dict(analysis=label, eICU_mean_SOFA=E[sofa].mean(), eICU_renal_ge3=((E.sofa_dynamic_renal_max.astype(float) >= 3) if sofa == 'sofa' else (E[sofa] - E.sofa + E.sofa_dynamic_renal_max.astype(float) >= 3)).mean(),
               AUC0=m0['AUC'], AUC1=m1['AUC'], dAUC=m1['AUC'] - m0['AUC'], Brier0=m0['Brier'], Brier1=m1['Brier'], dBrier=m0['Brier'] - m1['Brier'],
               OE0=m0['OE'], OE1=m1['OE'], Int0=m0['Intercept'], Int1=m1['Intercept'], Slope0=m0['Slope'], Slope1=m1['Slope'],
               DeLongP=delong(ye, 1 / (1 + np.exp(-lp1)), 1 / (1 + np.exp(-lp0))))
    for k, v in B.items(): row[k + '_lo'], row[k + '_hi'] = np.percentile(v, [2.5, 97.5])
    out.append(row)
pd.DataFrame(out).to_csv(OD + '/renal_sofa.csv', index=False)
print(pd.DataFrame(out).T.to_string())
