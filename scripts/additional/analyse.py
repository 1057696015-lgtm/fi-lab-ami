"""Additional analyses (NumPy only): cohort associations and splines (Table S5, Figure S3), FI-Lab values (Figure S1),
external sensitivity analyses (Figure 3, Table S6), severity-score comparisons and overlap
(Figure 3, Table S10), component observation/abnormality (Table S4, Figure S2).
Usage: python -I analyse.py OUT/data/cohorts.pkl OUT/results/additional
"""
import sys, json, math, numpy as np, pandas as pd
PK, OD = sys.argv[1:3]
res = pd.read_pickle(PK)
B = 1000

def logit(X, y, w=None):
    b = np.zeros(X.shape[1])
    for _ in range(100):
        eta = np.clip(X @ b, -30, 30); p = 1 / (1 + np.exp(-eta))
        g = X.T @ (y - p); H = (X * (p * (1 - p))[:, None]).T @ X
        s = np.linalg.solve(H, g); b += s
        if np.abs(s).max() < 1e-10: break
    return b

def robust(X, y, b):
    p = 1 / (1 + np.exp(-X @ b)); H = (X * (p * (1 - p))[:, None]).T @ X; Hi = np.linalg.inv(H)
    u = X * (y - p)[:, None]; n, k = X.shape
    return Hi @ (u.T @ u) @ Hi * n / (n - k)

def auc(y, p):
    r = pd.Series(p).rank().values; n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n0 * n1)

def calib(y, lp):
    a = 0.
    for _ in range(100):
        q = 1 / (1 + np.exp(-(lp + a))); st = (y - q).sum() / (q * (1 - q)).sum(); a += st
        if abs(st) < 1e-12: break
    s = logit(np.column_stack([np.ones(len(y)), lp]), y)[1]
    return a, s

def metrics(y, lp):
    p = 1 / (1 + np.exp(-lp)); a, s = calib(y, lp)
    return dict(AUC=auc(y, p), Brier=np.mean((p - y) ** 2), OE=y.sum() / p.sum(), CITL=a, Slope=s)

def pnorm2(z): return math.erfc(abs(z) / math.sqrt(2))

def ci(a): return np.percentile(a, [2.5, 97.5])

def design(d, fi=None, reduced=False, sofa='sofa', cci='charlson', drop_sofa=False):
    cols = [np.ones(len(d)), d.age, d.sexM, d['eth_Black'], d['eth_Asian'], d['eth_Other/Unknown'], d.bmi30, d.bmiM]
    if not reduced:
        cols += [d[c] * 1. for c in ['hypertension', 'diabetes', 'congestive_heart_failure', 'severe_liver_disease', 'renal_disease', 'chronic_pulmonary_disease']] + [d[cci]]
    if not drop_sofa: cols.append(d[sofa])
    if fi is not None: cols.append(fi)
    return np.column_stack(cols).astype(float)

def cohort(co, nmin=18, nmax=99):
    d = res[co]['d']; return d[(d.ncomp >= nmin) & (d.ncomp <= nmax) & d.sofa.notna() & (d.fol > 1)].copy()

M, E = cohort('M'), cohort('E')
assert len(M) == 4723 and len(E) == 4223
out = {'n': {'M': len(M), 'E': len(E), 'M_ev': int(M.death30.sum()), 'E_ev': int(E.death30.sum())}}

def fi_from(co, d, drop=()):
    F = res[co]['F'].loc[d.index].drop(columns=list(drop)); return (F.sum(1) / F.notna().sum(1)).values

# ---------- external frozen comparison helper ----------
def external(name, Md, Ed, kw0={}, kw1=None, fiM=None, fiE=None, seed=20261029):
    kw1 = kw0 if kw1 is None else kw1
    fiM = Md.fi_any.values * 10 if fiM is None else fiM * 10
    fiE = Ed.fi_any.values * 10 if fiE is None else fiE * 10
    ym, ye = Md.death30.values, Ed.death30.values
    b0 = logit(design(Md, **kw0), ym); b1 = logit(design(Md, fi=fiM, **kw1), ym)
    lp0 = design(Ed, **kw0) @ b0; lp1 = design(Ed, fi=fiE, **kw1) @ b1
    m0, m1 = metrics(ye, lp0), metrics(ye, lp1)
    rng = np.random.default_rng(seed); bs = []
    for _ in range(B):
        i = rng.integers(0, len(ye), len(ye)); yy = ye[i]
        if yy.sum() in (0, len(yy)): continue
        a0, a1 = metrics(yy, lp0[i]), metrics(yy, lp1[i])
        bs.append([a0['AUC'], a1['AUC'], a1['AUC'] - a0['AUC'], a0['Brier'] - a1['Brier'], a0['OE'], a1['OE'], a0['Slope'], a1['Slope']])
    bs = np.array(bs)
    r = dict(analysis=name, nM=len(Md), evM=int(ym.sum()), nE=len(Ed), evE=int(ye.sum()),
             auc0=m0['AUC'], auc1=m1['AUC'], dauc=m1['AUC'] - m0['AUC'], dbrier=m0['Brier'] - m1['Brier'],
             oe0=m0['OE'], oe1=m1['OE'], slope0=m0['Slope'], slope1=m1['Slope'], citl0=m0['CITL'], citl1=m1['CITL'], brier0=m0['Brier'], brier1=m1['Brier'])
    for j, k in enumerate(['auc0', 'auc1', 'dauc', 'dbrier', 'oe0', 'oe1', 'slope0', 'slope1']):
        r[k + '_lo'], r[k + '_hi'] = ci(bs[:, j])
    r['nboot'] = len(bs)
    return r, lp0, lp1

ext = []
r, lp0_main, lp1_main = external('Main analysis', M, E); ext.append(r)
ext.append(external('Creatinine-based SOFA renal subscore', M, E, kw0=dict(sofa='sofa_cr'))[0])
ext.append(external('Comorbidities and CCI removed', M, E, kw0=dict(reduced=True))[0])
ext.append(external('Unharmonized eICU CCI', M, E.assign(charlson=E.cci_raw))[0])
for lab, drop in [('FI-Lab without troponin T', ['Troponin_T']), ('FI-Lab without LDH', ['LDH']),
                  ('FI-Lab without PaO2, pH, lactate and troponin T', ['PaO2', 'pH', 'Lactate', 'Troponin_T']),
                  ('FI-Lab without SOFA-related items', ['Platelets', 'Bilirubin', 'Creatinine', 'PaO2', 'SBP', 'DBP'])]:
    ext.append(external(lab, M, E, fiM=fi_from('M', M, drop), fiE=fi_from('E', E, drop))[0])
# below availability threshold (12-17 components), models frozen from main development cohort
pd.DataFrame(ext).to_csv(OD + '/external_sensitivity.csv', index=False)

# ---------- within-database comparison with APACHE IV (eICU) and APS III (MIMIC-IV) ----------
def within(name, d, score, seed=11):
    d = d[d.sev.notna()]; y = d.death30.values; s = d.sev.values.astype(float); f = d.fi_any.values * 10
    X0 = np.column_stack([np.ones(len(d)), s]); X1 = np.column_stack([X0, f])
    def fitauc(Xa, Xb, yy, Ta, Tb, yt):
        ba, bb = logit(Xa, yy), logit(Xb, yy)
        return auc(yt, Ta @ ba), auc(yt, Tb @ bb)
    a0, a1 = fitauc(X0, X1, y, X0, X1, y)
    rng = np.random.default_rng(seed); bs, opt = [], []
    for _ in range(B):
        i = rng.integers(0, len(y), len(y))
        b0, b1 = logit(X0[i], y[i]), logit(X1[i], y[i])
        ab0, ab1 = auc(y[i], X0[i] @ b0), auc(y[i], X1[i] @ b1)
        ao0, ao1 = auc(y, X0 @ b0), auc(y, X1 @ b1)
        bs.append([ab0, ab1, ab1 - ab0]); opt.append([ab0 - ao0, ab1 - ao1, (ab1 - ab0) - (ao1 - ao0)])
    bs, opt = np.array(bs), np.array(opt); o = opt.mean(0)
    # frozen clinical models in same patients (eICU only, from main development)
    r = dict(analysis=name, score=score, n=len(d), ev=int(y.sum()), auc0=a0, auc1=a1, dauc=a1 - a0,
             auc0_lo=ci(bs[:, 0])[0], auc0_hi=ci(bs[:, 0])[1], auc1_lo=ci(bs[:, 1])[0], auc1_hi=ci(bs[:, 1])[1],
             dauc_lo=ci(bs[:, 2])[0], dauc_hi=ci(bs[:, 2])[1], auc0_oc=a0 - o[0], auc1_oc=a1 - o[1], dauc_oc=(a1 - a0) - o[2])
    return r, d
wa, Ea = within('eICU, APACHE IV', E, 'APACHE IV')
wb, Mb = within('MIMIC-IV, APS III', M, 'APS III', seed=12)
# frozen clinical models in APACHE subset
ix = E.index.get_indexer(Ea.index)
wa['frozen_auc0'] = auc(Ea.death30.values, lp0_main[ix]); wa['frozen_auc1'] = auc(Ea.death30.values, lp1_main[ix])
# APACHE IV + full clinical + SOFA, refitted within eICU, +FI (apparent & optimism)
def within_full(d, seed=13):
    d = d[d.sev.notna()]; y = d.death30.values
    X0 = np.column_stack([design(d), d.sev.values]); X1 = np.column_stack([X0, d.fi_any.values * 10])
    b0, b1 = logit(X0, y), logit(X1, y); a0, a1 = auc(y, X0 @ b0), auc(y, X1 @ b1)
    rng = np.random.default_rng(seed); bs, opt = [], []
    for _ in range(B):
        i = rng.integers(0, len(y), len(y)); c0, c1 = logit(X0[i], y[i]), logit(X1[i], y[i])
        ab0, ab1 = auc(y[i], X0[i] @ c0), auc(y[i], X1[i] @ c1); bs.append(ab1 - ab0)
        opt.append((ab1 - ab0) - (auc(y, X1 @ c1) - auc(y, X0 @ c0)))
    return dict(auc0=a0, auc1=a1, dauc=a1 - a0, dauc_lo=ci(bs)[0], dauc_hi=ci(bs)[1], dauc_oc=a1 - a0 - np.mean(opt))
wa['full'] = within_full(E)
json.dump({'apache': wa, 'aps': wb}, open(OD + '/within_severity.json', 'w'), indent=1, default=float)

# ---------- associations ----------
def assoc(d):
    y = d.death30.values; f = d.fi_any.values * 10; o = {}
    for lab, X in [('unadj', np.column_stack([np.ones(len(d)), f])),
                   ('demo', np.column_stack([np.ones(len(d)), d.age, d.sexM, d['eth_Black'], d['eth_Asian'], d['eth_Other/Unknown'], d.bmi30, d.bmiM, f]).astype(float)),
                   ('nosofa', design(d, fi=f, drop_sofa=True)),
                   ('full', design(d, fi=f))]:
        b = logit(X, y); V = robust(X, y, b); se = math.sqrt(V[-1, -1])
        o[lab] = dict(b=b[-1], se=se, OR=math.exp(b[-1]), lo=math.exp(b[-1] - 1.959964 * se), hi=math.exp(b[-1] + 1.959964 * se), p=pnorm2(b[-1] / se))
    return o
A = {'M': assoc(M), 'E': assoc(E)}
for lab in A['M']:
    bm, be = A['M'][lab], A['E'][lab]; z = (be['b'] - bm['b']) / math.sqrt(be['se'] ** 2 + bm['se'] ** 2); sd = math.sqrt(be['se'] ** 2 + bm['se'] ** 2)
    A.setdefault('ratio', {})[lab] = dict(R=math.exp(be['b'] - bm['b']), lo=math.exp(be['b'] - bm['b'] - 1.959964 * sd), hi=math.exp(be['b'] - bm['b'] + 1.959964 * sd), p=pnorm2(z))

# spline (natural cubic, 4 knots at min, 33.3, 66.7 pct, max == R ns(df=3)); reference 0.30
def rcs_basis(x, k):
    def tp(u): return np.maximum(u, 0) ** 3
    K = len(k); cols = [x]
    for j in range(K - 2):
        cols.append(tp(x - k[j]) - tp(x - k[K - 2]) * (k[K - 1] - k[j]) / (k[K - 1] - k[K - 2]) + tp(x - k[K - 1]) * (k[K - 2] - k[j]) / (k[K - 1] - k[K - 2]))
    return np.column_stack(cols)
spl = {}
for co, d in [('M', M), ('E', E)]:
    x = d.fi_any.values; y = d.death30.values
    k = np.r_[x.min(), np.quantile(x, [1 / 3, 2 / 3]), x.max()]
    C = design(d); X = np.column_stack([C, rcs_basis(x, k)]); b = logit(X, y); V = robust(X, y, b)
    nb = 3; bs, Vs = b[-nb:], V[-nb:, -nb:]
    nl = bs[1:]; W = float(nl @ np.linalg.solve(Vs[1:, 1:], nl)); pnl = math.exp(-W / 2)
    lo_, hi_ = np.quantile(x, [.01, .99]); grid = np.linspace(max(lo_, .10), min(hi_, .90), 120)
    G = rcs_basis(grid, k) - rcs_basis(np.array([.30]), k); est = G @ bs; se = np.sqrt(np.einsum('ij,jk,ik->i', G, Vs, G))
    spl[co] = dict(grid=grid.tolist(), lor=est.tolist(), se=se.tolist(), p_nonlin=pnl, knots=k.tolist())
json.dump({'assoc': A, 'spline': spl}, open(OD + '/associations.json', 'w'), indent=1, default=float)

# ---------- overlap ----------
def spearman(a, b):
    a, b = pd.Series(a).rank().values, pd.Series(b).rank().values; return np.corrcoef(a, b)[0, 1]
ov = {}
for co, d in [('M', M), ('E', E)]:
    dd = d[d.sev.notna()]
    r = dict(n=len(dd), sev=spearman(dd.fi_any, dd.sev), sofa=spearman(d.fi_any, d.sofa), cci=spearman(d.fi_any, d.charlson), age=spearman(d.fi_any, d.age))
    X = np.column_stack([np.ones(len(dd)), dd.sev, dd.sofa, dd.charlson, dd.age]).astype(float); bb = np.linalg.lstsq(X, dd.fi_any.values, rcond=None)[0]
    res_ = dd.fi_any.values - X @ bb; r['R2'] = 1 - res_.var() / dd.fi_any.values.var()
    r['ncomp_med'] = float(d.ncomp.median()); r['ncomp_q'] = d.ncomp.quantile([.25, .75]).tolist()
    r['sev_med'] = float(dd.sev.median()); r['sev_q'] = dd.sev.quantile([.25, .75]).tolist(); r['sev_avail'] = len(dd)
    ov[co] = r
json.dump(ov, open(OD + '/overlap.json', 'w'), indent=1, default=float)

# ---------- component prevalence and raw medians ----------
rows = []
for co, d in [('M', M), ('E', E)]:
    F = res[co]['F'].loc[d.index]; V = res[co]['V'].loc[d.index]
    for it in F.columns:
        rows.append(dict(cohort=co, item=it, observed=F[it].notna().mean(), abnormal=F[it].mean(), median=V[it].median(), q1=V[it].quantile(.25), q3=V[it].quantile(.75)))
pd.DataFrame(rows).to_csv(OD + '/components.csv', index=False)

# save predictions for figures
pd.DataFrame(dict(y=E.death30.values, lp0=lp0_main, lp1=lp1_main, fi=E.fi_any.values)).to_csv(OD + '/eicu_pred_main.csv', index=False)
pd.DataFrame(dict(co=['M'] * len(M) + ['E'] * len(E), fi=np.r_[M.fi_any.values, E.fi_any.values])).to_csv(OD + '/fi_values.csv', index=False)
print(pd.DataFrame(ext)[['analysis', 'nE', 'evE', 'auc0', 'auc1', 'dauc', 'dauc_lo', 'dauc_hi', 'oe1', 'slope1']].round(3).to_string())
print(json.dumps({'apache': wa, 'aps': wb}, default=float, indent=0)[:1500])
print(json.dumps(A, default=float)[:2000]); print({k: v['p_nonlin'] for k, v in spl.items()}); print(ov)
