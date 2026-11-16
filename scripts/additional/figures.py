"""Figures 1-3 and Supplementary Figures S1-S4.
Usage: python -I figures.py ADDITIONAL_RESULTS MAIN_RESULTS COHORTS_PKL FIGDIR COHORT_FLOW_JSON
MAIN_RESULTS must contain external_validation.csv, internal_validation.csv and decision_curves*.csv from the R main analysis."""
import sys, json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
RD, REF, PK, FD = sys.argv[1:5]
FLOW = sys.argv[5]
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.titleweight': 'bold', 'axes.titlesize': 12, 'axes.titlelocation': 'left', 'savefig.dpi': 600})
CM, CE, GR = '#2B6A94', '#B8642E', '#5A6470'

def lab(ax, s, x=-0.12, y=1.06): ax.text(x, y, s, transform=ax.transAxes, fontsize=14, fontweight='bold', va='bottom')

# ---------------- Figure 1: (a) study design, (b) cohort flow ----------------
FL = {co: dict(v) for co, v in json.load(open(FLOW)).items()}
fig = plt.figure(figsize=(7.2, 8.8))
ax = fig.add_axes([0, 0.705, 1, 0.295]); ax.axis('off'); ax.set_xlim(0, 100); ax.set_ylim(0, 40)
ax.text(1.5, 38.5, 'a', fontsize=14, fontweight='bold', va='top')
nM, eM, nE, eE = FL['M']['Analysis cohort'], FL['M']['Analysis cohort deaths'], FL['E']['Analysis cohort'], FL['E']['Analysis cohort deaths']
for x0, c, fill, l1, l2, l3, l4 in [(4, CM, '#EAF1F6', 'MIMIC-IV v3.1', 'Single center, 2008–2022', 'Model development', f'n = {nM:,}; {eM} deaths'),
                                    (66, CE, '#F8EEE6', 'eICU-CRD', 'Multicenter, 2014–2015', 'External validation', f'n = {nE:,}; {eE} deaths')]:
    ax.add_patch(FancyBboxPatch((x0, 21), 30, 17, boxstyle='round,pad=0,rounding_size=1.2', fc=fill, ec=c, lw=1.2))
    ax.text(x0 + 15, 34.6, l1, ha='center', va='center', fontsize=10.5, fontweight='bold', color=c)
    ax.text(x0 + 15, 31.0, l2, ha='center', va='center', fontsize=9)
    ax.text(x0 + 15, 27.6, l3, ha='center', va='center', fontsize=9, fontweight='bold')
    ax.text(x0 + 15, 24.2, l4, ha='center', va='center', fontsize=9)
ax.annotate('', xy=(65.4, 29.5), xytext=(34.6, 29.5), arrowprops=dict(arrowstyle='-|>', color='#333', lw=1.4, mutation_scale=14))
ax.text(50, 33.6, 'Clinical + SOFA with and without FI-Lab', ha='center', va='center', fontsize=8.2, fontweight='bold')
ax.text(50, 31.2, 'coefficients frozen', ha='center', va='center', fontsize=8.6)
ax.text(50, 27.4, 'Same 27-item FI-Lab', ha='center', va='center', fontsize=8.6)
ax.text(50, 25.2, 'and thresholds', ha='center', va='center', fontsize=8.6)
# timeline: -6 h ... 24 h (linear), break, day 30
X = lambda h: 10 + (h + 6) * (52 / 30); yT = 9
ax.add_patch(plt.Rectangle((X(-6), yT - 1.6), X(24) - X(-6), 3.2, fc='#E9EEF2', ec='none'))
ax.plot([X(-6), X(24)], [yT, yT], color='#333', lw=1.3); ax.plot([X(24), 84], [yT, yT], color='#333', lw=1.3)
ax.plot([84.6, 85.6], [yT - 1.2, yT + 1.2], color='#333', lw=1); ax.plot([85.6, 86.6], [yT - 1.2, yT + 1.2], color='#333', lw=1)
ax.plot([86.6, 93], [yT, yT], color='#333', lw=1.3)
for xx, lab_ in [(X(-6), '−6 h'), (X(0), 'ICU admission'), (X(24), '24 h'), (93, 'Day 30')]:
    ax.plot([xx, xx], [yT - 1.6, yT + 1.6], color='#333', lw=1.1); ax.text(xx, yT - 3.6, lab_, ha='center', va='center', fontsize=8.6)
ax.text((X(-6) + X(24)) / 2, yT + 4.2, 'Predictor window: FI-Lab (first values), SOFA', ha='center', va='center', fontsize=8.8)
ax.plot(X(24), yT, marker='v', color=CE, ms=9, transform=ax.transData); ax.plot([X(24), X(24)], [yT + 1.6, yT + 6.6], color=CE, lw=1)
ax.text(X(24), yT + 8.6, 'Prediction point: end of first ICU day (alive and hospitalized)', ha='center', va='center', fontsize=8.8, color=CE, fontweight='bold')
ax.text((X(24) + 93) / 2 + 1, yT + 4.2, 'Outcome: in-hospital death by day 30', ha='center', va='center', fontsize=8.8)
ax = fig.add_axes([0, 0, 1, 0.70]); ax.axis('off'); ax.set_xlim(0, 100); ax.set_ylim(0, 100)
ax.text(1.5, 99, 'b', fontsize=14, fontweight='bold', va='top')
cols = {'M': (4, 46, CM, '#EAF1F6', 'MIMIC-IV development'), 'E': (54, 96, CE, '#F8EEE6', 'eICU external validation')}
def steps_for(co):
    f = FL[co]; c0 = f['Candidate AMI ICU records']; c1 = f['Adult, first ICU record per hospitalization']; c2 = f['One admission per patient']
    c3 = f['Valid hospital outcome and timing']; c4 = f['Complete covariates']; c5 = f['FI-Lab >=18 components']; c6 = f['FI-Lab >=18 components and SOFA available']
    dx, ox, n, ev = f['Excluded within 24 hours: deaths'], f['Excluded within 24 hours: other discharges'], f['Analysis cohort'], f['Analysis cohort deaths']
    rule = 'earliest' if co == 'M' else 'random'
    return [('Candidate AMI ICU records', f'n = {c0:,}'), ('Adult, first ICU record per hospitalization', f'n = {c1:,}; excluded = {c0 - c1:,}'),
            (f'One admission per patient ({rule})', f'n = {c2:,}; repeat admissions excluded = {c1 - c2:,}'), ('Valid hospital outcome and timing', f'n = {c3:,}; excluded = {c2 - c3:,}'),
            ('FI-Lab ≥18 components and SOFA available', f'n = {c6:,}; <18 components = {c4 - c5:,}; no SOFA = {c5 - c6:,}'),
            ('Alive and hospitalized beyond 24 hours', f'Excluded: {dx} deaths, {ox} discharges'),
            ('Development cohort' if co == 'M' else 'External validation cohort', f'n = {n:,}; {ev} deaths ({100 * ev / n:.1f}%)')]
H, gap, top = 9.6, 3.2, 91
for co, (x0, x1, c, fill, title) in cols.items():
    ax.text((x0 + x1) / 2, 95.5, title, ha='center', va='center', fontsize=11.5, fontweight='bold', color=c)
    st = steps_for(co)
    for i, (a_, b_) in enumerate(st):
        y = top - i * (H + gap) - H; last = i == len(st) - 1
        ax.add_patch(FancyBboxPatch((x0, y), x1 - x0, H, boxstyle='round,pad=0,rounding_size=1.0', fc=fill if last else 'white', ec=c, lw=1.1))
        ax.text((x0 + x1) / 2, y + H * 0.66, a_, ha='center', va='center', fontsize=9, fontweight='bold' if last else 'normal')
        ax.text((x0 + x1) / 2, y + H * 0.30, b_, ha='center', va='center', fontsize=8.7)
        if not last: ax.annotate('', xy=((x0 + x1) / 2, y - gap + 0.3), xytext=((x0 + x1) / 2, y - 0.3), arrowprops=dict(arrowstyle='->', color='#777', lw=1))
fig.savefig(FD + '/Figure1_design_flow.png'); plt.close(fig)

# ---------------- Figure 2: frozen-model external validation ----------------
pr = pd.read_csv(RD + '/eicu_pred_main.csv'); y = pr.y.values.astype(int)
ext = pd.read_csv(REF + '/external_validation.csv').set_index('metric'); inn = pd.read_csv(REF + '/internal_validation.csv').set_index('metric')
ren = pd.read_csv(RD + '/renal_sofa.csv').set_index('analysis').loc['Creatinine-based renal SOFA']
fig = plt.figure(figsize=(7.2, 6.3))
axa = fig.add_axes([0.09, 0.42, 0.37, 0.50]); axb = fig.add_axes([0.60, 0.42, 0.37, 0.50]); axc = fig.add_axes([0.40, 0.08, 0.36, 0.22])
for k, c, ls, nm in [('lp0', GR, '-', 'Clinical + SOFA'), ('lp1', CE, '--', 'Clinical + SOFA + FI-Lab')]:
    s = pr[k].values; o = np.argsort(-s, kind='mergesort'); ys = y[o]
    tpr = np.r_[0, np.cumsum(ys) / ys.sum()]; fpr = np.r_[0, np.cumsum(1 - ys) / (1 - ys).sum()]
    a_ = ext.loc['Without.AUC' if k == 'lp0' else 'With.AUC', 'estimate']
    axa.plot(fpr, tpr, color=c, ls=ls, lw=1.6, label=f'{nm} (AUC {a_:.3f})')
    p_ = 1 / (1 + np.exp(-s)); grp = np.minimum(10, np.ceil(pd.Series(p_).rank(method='first').values / len(p_) * 10)).astype(int)
    P, Ob, Lo, Hi = [], [], [], []
    for g_ in range(1, 11):
        ii = grp == g_; n = ii.sum(); ob = y[ii].mean(); z = 1.959964; cen = (ob + z * z / (2 * n)) / (1 + z * z / n)
        hw = z * np.sqrt(ob * (1 - ob) / n + z * z / (4 * n * n)) / (1 + z * z / n); P.append(p_[ii].mean()); Ob.append(ob); Lo.append(cen - hw); Hi.append(cen + hw)
    P, Ob = np.array(P), np.array(Ob)
    axb.errorbar(P, Ob, yerr=[Ob - np.array(Lo), np.array(Hi) - Ob], color=c, ls=ls, lw=1.3, marker='o', ms=3.5, capsize=2, elinewidth=0.9)
axa.plot([0, 1], [0, 1], color='#BBB', lw=0.8); axa.set_xlim(0, 1); axa.set_ylim(0, 1); axa.set_aspect('equal')
axa.set_xlabel('1 − Specificity'); axa.set_ylabel('Sensitivity'); axa.set_title('External discrimination'); axa.legend(frameon=False, loc='lower right', fontsize=8)
axb.plot([0, 0.6], [0, 0.6], color='#BBB', lw=0.8); axb.set_xlim(0, 0.6); axb.set_ylim(0, 0.6); axb.set_aspect('equal')
axb.set_xticks(np.arange(0, 0.61, 0.1)); axb.set_yticks(np.arange(0, 0.61, 0.1))
axb.set_xlabel('Mean predicted risk'); axb.set_ylabel('Observed mortality'); axb.set_title('External calibration')
lab(axa, 'a', -0.2, 1.04); lab(axb, 'b', -0.2, 1.04)
rows_c = [('MIMIC-IV internal (optimism-corrected)', inn.loc['Delta_AUC', 'optimism_corrected'], None, None, CM),
          ('eICU external, exported SOFA', ext.loc['Delta_AUC', 'estimate'], ext.loc['Delta_AUC', 'lower'], ext.loc['Delta_AUC', 'upper'], '#1F1F1F'),
          ('eICU external, creatinine-based renal SOFA', ren.dAUC, ren.dAUC_lo, ren.dAUC_hi, CE)]
import matplotlib.transforms as mt
BTc = mt.blended_transform_factory(fig.transFigure, axc.transData)
for i, (nm, est, lo, hi, c) in enumerate(rows_c):
    yy = len(rows_c) - 1 - i
    if lo is None:
        axc.plot(est, yy, 'o', mfc='white', mec=c, mew=1.6, ms=7); txt = f'{est:.3f}'
    else:
        axc.plot([lo, hi], [yy, yy], color=c, lw=1.6); axc.plot(est, yy, 'o', color=c, ms=6); txt = f'{est:.3f} ({lo:.3f} to {hi:.3f})'
    axc.text(0.02, yy, nm, transform=BTc, fontsize=8.8, va='center', ha='left', clip_on=False)
    axc.text(0.78, yy, txt, transform=BTc, fontsize=8.8, va='center', ha='left', clip_on=False)
axc.axvline(0, color='#888', ls='--', lw=0.9); axc.set_xlim(-0.003, 0.05); axc.set_ylim(-0.6, len(rows_c) - 0.4)
axc.set_yticks([]); axc.spines['left'].set_visible(False); axc.set_xticks([0, 0.01, 0.02, 0.03, 0.04, 0.05]); axc.set_xlabel('Increase in AUC with FI-Lab')
fig.text(0.02, 0.33, 'c', fontsize=14, fontweight='bold', va='bottom'); fig.text(0.06, 0.33, 'Internal and external increments', fontsize=12, fontweight='bold', va='bottom')
fig.savefig(FD + '/Figure2_external_validation.png'); plt.close(fig)

# ---------------- Figure S1: distribution; Figure S3: spline ----------------
fv = pd.read_csv(RD + '/fi_values.csv'); sp = json.load(open(RD + '/associations.json'))['spline']
fig, ax = plt.subplots(figsize=(7.2, 3.4)); fig.subplots_adjust(left=0.09, right=0.98, bottom=0.16, top=0.97)
g = np.linspace(0, 1, 400)
for co, c, nm in [('M', CM, 'MIMIC-IV development'), ('E', CE, 'eICU validation')]:
    x = fv[fv.co == co].fi.values; bw = 1.06 * x.std() * len(x) ** -0.2
    den = np.exp(-0.5 * ((g[:, None] - x[None]) / bw) ** 2).sum(1) / (len(x) * bw * np.sqrt(2 * np.pi))
    ax.plot(g, den, color=c, lw=1.8, label=nm); ax.fill_between(g, den, color=c, alpha=0.12)
for cut in (0.40, 0.53): ax.axvline(cut, color='#555', ls=':', lw=1)
ax.text(0.465, 2.85, 'MIMIC-IV tertile cut-points', ha='center', va='top', fontsize=8, color='#555')
ax.set_xlim(0, 1); ax.set_ylim(0, 3.0); ax.set_xlabel('Common 27-item FI-Lab'); ax.set_ylabel('Density'); ax.legend(frameon=False, loc='upper right', fontsize=9)
fig.savefig(FD + '/FigureS1_distribution.png'); plt.close(fig)
fig, ax = plt.subplots(figsize=(4.6, 3.9)); fig.subplots_adjust(left=0.16, right=0.97, bottom=0.14, top=0.97)
for co, c, ls, nm in [('M', CM, '-', 'MIMIC-IV'), ('E', CE, '--', 'eICU')]:
    s = sp[co]; gr = np.array(s['grid']); l = np.array(s['lor']); se = np.array(s['se'])
    p = s['p_nonlin']; pt = 'P <0.001' if p < 0.001 else f'P = {p:.2f}'
    ax.plot(gr, np.exp(l), color=c, ls=ls, lw=1.8, label=f'{nm} (nonlinearity {pt})'); ax.fill_between(gr, np.exp(l - 1.96 * se), np.exp(l + 1.96 * se), color=c, alpha=0.15, lw=0)
ax.set_yscale('log'); ax.axhline(1, color='#888', ls=':', lw=1); ax.set_xlabel('Common 27-item FI-Lab'); ax.set_ylabel('Adjusted odds ratio (reference 0.30)')
ax.set_yticks([0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30]); ax.set_yticklabels(['0.01', '0.03', '0.1', '0.3', '1', '3', '10', '30']); ax.minorticks_off()
ax.legend(frameon=False, loc='upper left', fontsize=8.5)
fig.savefig(FD + '/FigureS3_spline.png'); plt.close(fig)

# ---------------- Figure 3: forest (comparator, robustness, eICU regions) ----------------
ex = pd.read_csv(RD + '/external_sensitivity.csv').set_index('analysis'); ws = json.load(open(RD + '/within_severity.json'))
ref = pd.read_csv(REF + '/external_validation.csv').set_index('metric')
rg = pd.read_csv(RD + '/eicu_region_validation.csv').set_index('region')
rs = pd.read_csv(RD + '/renal_sofa.csv').set_index('analysis')
rows = [('Comparator', None)]
def add(name, est, lo, hi, n, ev, kind='E'): rows.append((name, (est, lo, hi, n, ev, kind)))
add('Clinical + exported SOFA (main; frozen)', ref.loc['Delta_AUC', 'estimate'], ref.loc['Delta_AUC', 'lower'], ref.loc['Delta_AUC', 'upper'], 4223, 451, 'main')
r = rs.loc['Creatinine-based renal SOFA']; add('Clinical + creatinine-based SOFA (frozen)', r.dAUC, r.dAUC_lo, r.dAUC_hi, 4223, 451)
b = ws['aps']; add('APS III (MIMIC-IV, refitted)', b['dauc'], b['dauc_lo'], b['dauc_hi'], b['n'], b['ev'], 'M')
a = ws['apache']; add('APACHE IV (eICU, refitted)', a['dauc'], a['dauc_lo'], a['dauc_hi'], a['n'], a['ev'])
rows.append(('Robustness of the main comparison (frozen)', None))
for k, nm in [('Comorbidities and CCI removed', 'Comorbidities and CCI removed'), ('Unharmonized eICU CCI', 'Unharmonized eICU CCI'),
              ('FI-Lab without troponin T', 'FI-Lab without troponin T'), ('FI-Lab without LDH', 'FI-Lab without LDH'),
              ('FI-Lab without PaO2, pH, lactate and troponin T', 'FI-Lab without PaO$_2$, pH, lactate, troponin T'),
              ('FI-Lab without SOFA-related items', 'FI-Lab without SOFA-overlapping items')]:
    r = ex.loc[k]; add(nm, r.dauc, r.dauc_lo, r.dauc_hi, int(r.nE), int(r.evE))
rows.append(('eICU census region (frozen main models)', None))
for k in ['Midwest', 'Northeast', 'South', 'West', 'Missing']:
    r = rg.loc[k]; add('Region not recorded' if k == 'Missing' else k, r.dAUC, r.dAUC_lo, r.dAUC_hi, int(r.n), int(r.deaths))
fig = plt.figure(figsize=(7.2, 7.3)); ax = fig.add_axes([0.43, 0.075, 0.36, 0.905])
import matplotlib.transforms as mt
BT = mt.blended_transform_factory(fig.transFigure, ax.transData)
y = len(rows)
for nm, v in rows:
    y -= 1
    if v is None:
        ax.text(0.01, y, nm, fontsize=9.5, fontweight='bold', ha='left', va='center', transform=BT, clip_on=False); continue
    est, lo, hi, n, ev, kind = v; c = {'main': '#1F1F1F', 'M': CM}.get(kind, CE)
    ax.plot([max(lo, -0.03), min(hi, 0.12)], [y, y], color=c, lw=1.6); ax.plot(est, y, 'D' if kind == 'main' else 'o', color=c, ms=6 if kind == 'main' else 5)
    ax.text(0.025, y + 0.16, nm, fontsize=8.6, ha='left', va='center', transform=BT, clip_on=False)
    ax.text(0.025, y - 0.22, f'n = {n:,}; deaths = {ev}', fontsize=7.6, color='#555', ha='left', va='center', transform=BT, clip_on=False)
    ax.text(0.80, y, f'{est:.3f} ({lo:.3f} to {hi:.3f})'.replace('-', '−'), fontsize=8.6, va='center', transform=BT, clip_on=False)
ax.axvline(0, color='#888', ls='--', lw=0.9); ax.set_xlim(-0.03, 0.12); ax.set_ylim(-0.7, len(rows) - 0.4)
ax.set_yticks([]); ax.spines['left'].set_visible(False); ax.set_xlabel('Increase in AUC with FI-Lab (95% CI)')
ax.set_xticks([-0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.10]); ax.set_xticklabels(['−0.02', '0', '0.02', '0.04', '0.06', '0.08', '0.10'])
fig.savefig(FD + '/Figure3_forest.png'); plt.close(fig)

# ---------------- Figure S4: decision curves ----------------
import glob
dc = pd.read_csv(sorted(glob.glob(REF + '/decision_curves*.csv'))[0])
dc = dc[dc.scenario == 'survivor24'] if 'scenario' in dc else dc
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.5)); fig.subplots_adjust(left=0.09, right=0.98, bottom=0.15, top=0.86, wspace=0.34)
ax = axs[0]; t = dc.threshold
ax.plot(t, dc['all'], color='#AAA', lw=1.1, label='Treat all'); ax.axhline(0, color='#444', lw=0.9, label='Treat none')
ax.plot(t, dc.without, color=GR, lw=1.8, label='Clinical + SOFA'); ax.plot(t, dc['with'], color=CE, lw=1.8, ls='--', label='Clinical + SOFA + FI-Lab')
ax.set_ylim(-0.01, 0.085); ax.set_xlim(0.05, 0.5); ax.set_xlabel('Threshold probability'); ax.set_ylabel('Net benefit'); ax.set_title('Net benefit in eICU'); ax.legend(frameon=False, fontsize=8)
lab(ax, 'a')
ax = axs[1]; ax.fill_between(t, dc.lower, dc.upper, color=CE, alpha=0.18, lw=0); ax.plot(t, dc.difference, color=CE, lw=1.8)
ax.axhline(0, color='#888', ls='--', lw=0.9); ax.set_xlim(0.05, 0.5); ax.set_xlabel('Threshold probability'); ax.set_ylabel('Difference in net benefit')
ax.set_title('Net benefit added by FI-Lab'); lab(ax, 'b')
fig.savefig(FD + '/FigureS4_decision_curves.png'); plt.close(fig)

# ---------------- Figure S2: component abnormal prevalence ----------------
cp = pd.read_csv(RD + '/components.csv'); pv = cp.pivot(index='item', columns='cohort', values='abnormal'); ob = cp.pivot(index='item', columns='cohort', values='observed')
pv = pv.loc[ob.mean(1).sort_values(kind='mergesort').index]
fig, axs = plt.subplots(1, 2, figsize=(7.2, 6.4), sharey=True); fig.subplots_adjust(left=0.16, right=0.98, bottom=0.08, top=0.88, wspace=0.12)
yy = np.arange(len(pv)); nm = [i.replace('_', ' ').replace('Troponin T', 'Troponin T').replace('PaO2', 'PaO$_2$').replace('PaCO2', 'PaCO$_2$') for i in pv.index]
for ax, D, ttl in [(axs[0], pv, 'Abnormal among observed'), (axs[1], ob.loc[pv.index], 'Observed')]:
    ax.hlines(yy, D.min(1), D.max(1), color='#BBB', lw=1); ax.plot(D.M, yy, 'o', color=CM, ms=5, label='MIMIC-IV'); ax.plot(D.E, yy, 's', color=CE, ms=4.5, label='eICU')
    ax.set_xlim(0, 1.02); ax.set_title(ttl); ax.set_xlabel('Proportion of patients')
axs[0].set_yticks(yy); axs[0].set_yticklabels(nm, fontsize=8.5); h_, l_ = axs[0].get_legend_handles_labels(); fig.legend(h_, l_, frameon=False, loc='upper center', bbox_to_anchor=(0.57, 0.995), ncol=2, fontsize=9); lab(axs[0], 'a', -0.2); lab(axs[1], 'b', -0.05)
fig.savefig(FD + '/FigureS2_components.png'); plt.close(fig)
print('ok')
