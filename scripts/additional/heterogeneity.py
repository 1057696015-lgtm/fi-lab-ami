"""Additional analyses (2026-11-03, NumPy only):
(1) FI-Lab within SOFA categories: cohort-specific logistic model with all model covariates, SOFA
    (continuous) and SOFA-category indicators, with a separate FI-Lab slope per category (OR per 0.1);
    LR tests of FI-Lab x SOFA category (df 2) and FI-Lab x continuous SOFA (df 1).
(2) eICU region-stratified frozen validation: AUC without/with FI-Lab, paired bootstrap dAUC, O/E.
Usage: python -I heterogeneity.py OUT/data/cohorts.pkl OUTDIR"""
import sys, json, math, numpy as np, pandas as pd
PK, OD = sys.argv[1:3]; res = pd.read_pickle(PK)
def cohort(co):
    d = res[co]['d']; return d[(d.ncomp >= 18) & d.sofa.notna() & (d.fol > 1)].copy()
COV = ['hypertension', 'diabetes', 'congestive_heart_failure', 'severe_liver_disease', 'renal_disease', 'chronic_pulmonary_disease']
def base(d):
    return [np.ones(len(d)), d.age, d.sexM, d['eth_Black'], d['eth_Asian'], d['eth_Other/Unknown'], d.bmi30, d.bmiM] + [d[c] * 1. for c in COV] + [d.charlson, d.sofa]
def design(d, fi=True):
    cols = base(d) + ([d.fi_any.values * 10] if fi else []); return np.column_stack(cols).astype(float)
def logit(X, y):
    b = np.zeros(X.shape[1])
    for _ in range(100):
        p = 1 / (1 + np.exp(-X @ b)); H = (X * (p * (1 - p))[:, None]).T @ X; s = np.linalg.solve(H, X.T @ (y - p)); b += s
        if np.abs(s).max() < 1e-10: break
    p = 1 / (1 + np.exp(-X @ b)); H = (X * (p * (1 - p))[:, None]).T @ X
    ll = np.sum(y * np.log(p) + (1 - y) * np.log(1 - p))
    return b, np.linalg.inv(H), ll
def auc(y, p):
    r = pd.Series(p).rank().values; n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n0 * n1)
P = lambda X, b: 1 / (1 + np.exp(-X @ b))
out = {}
M, E = cohort('M'), cohort('E')
# (1)
def sofa_cat(s): return np.select([s <= 3, s <= 7], ['SOFA 0–3', 'SOFA 4–7'], 'SOFA ≥8')
rows = []
for co, d in [('MIMIC-IV', M), ('eICU', E)]:
    y = d.death30.values.astype(float)
    X = design(d); b, V, ll1 = logit(X, y)
    Xi = np.column_stack([X, d.sofa.values * d.fi_any.values * 10]); bi, Vi, lli = logit(Xi, y)
    lr = 2 * (lli - ll1); p_int = math.erfc(math.sqrt(lr / 2))
    out[co] = dict(interaction_LR=lr, interaction_P=p_int, interaction_coef=bi[-1], interaction_se=math.sqrt(Vi[-1, -1]))
    sc = sofa_cat(d.sofa.values); cats = ['SOFA 0–3', 'SOFA 4–7', 'SOFA ≥8']
    # one model, all covariates + SOFA (continuous) + SOFA-category indicators, FI-Lab slope per category
    ind = [(sc == c).astype(float) for c in cats]
    Xc = np.column_stack(base(d) + ind[1:] + [d.fi_any.values * 10 * g for g in ind]).astype(float)
    Xr = np.column_stack(base(d) + ind[1:] + [d.fi_any.values * 10]).astype(float)
    bc, Vc, llc = logit(Xc, y); _, _, llr = logit(Xr, y)
    lr2 = 2 * (llc - llr); out[co]['category_interaction_LR_df2'] = lr2; out[co]['category_interaction_P'] = math.exp(-lr2 / 2)
    for j, c in enumerate(cats):
        k = bc[len(bc) - 3 + j]; se = math.sqrt(Vc[len(bc) - 3 + j, len(bc) - 3 + j]); m = sc == c
        rows.append(dict(cohort=co, sofa=c, n=int(m.sum()), deaths=int(y[m].sum()), OR=math.exp(k), lo=math.exp(k - 1.959964 * se), hi=math.exp(k + 1.959964 * se),
                         P=math.erfc(abs(k / se) / math.sqrt(2))))
pd.DataFrame(rows).to_csv(OD + '/within_sofa_or.csv', index=False)
# (2) region-stratified frozen validation in eICU
ym = M.death30.values.astype(float); b0 = logit(design(M, False), ym)[0]; b1 = logit(design(M), ym)[0]
E['p0'] = P(design(E, False), b0); E['p1'] = P(design(E), b1)
reg = E.region.fillna('Missing').values
rr = []
rng = np.random.default_rng(20261103)
for g in ['Midwest', 'Northeast', 'South', 'West', 'Missing']:
    d = E[reg == g]
    if len(d) == 0: continue
    y = d.death30.values.astype(float); p0, p1 = d.p0.values, d.p1.values
    a0, a1 = auc(y, p0), auc(y, p1); bs = []
    for _ in range(1000):
        i = rng.integers(0, len(y), len(y))
        if y[i].sum() in (0, len(y)): continue
        bs.append(auc(y[i], p1[i]) - auc(y[i], p0[i]))
    lo, hi = np.percentile(bs, [2.5, 97.5])
    rr.append(dict(region=g, n=len(d), deaths=int(y.sum()), AUC0=a0, AUC1=a1, dAUC=a1 - a0, dAUC_lo=lo, dAUC_hi=hi, OE0=y.sum() / p0.sum(), OE1=y.sum() / p1.sum()))
pd.DataFrame(rr).to_csv(OD + '/eicu_region_validation.csv', index=False)
# overall check
y = E.death30.values.astype(float); out['check_overall'] = dict(AUC0=auc(y, E.p0.values), AUC1=auc(y, E.p1.values))
json.dump(out, open(OD + '/heterogeneity.json', 'w'), indent=1, ensure_ascii=False)
print(json.dumps(out, indent=1, ensure_ascii=False)); print(pd.DataFrame(rows).round(3).to_string()); print(pd.DataFrame(rr).round(3).to_string())
