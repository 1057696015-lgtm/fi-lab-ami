"""Figure 4: (a, b) observed 30-day mortality by SOFA category and FI-Lab tertile (cut-points from MIMIC-IV);
(c) adjusted OR per 0.1 FI-Lab within each SOFA category (heterogeneity.py).
Usage: python -I figure5.py OUT/results/additional OUT/figures"""
import sys, json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
RD, FD = sys.argv[1:3]
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 10, 'savefig.dpi': 600})
g = pd.read_csv(RD + '/risk_grid.csv'); c = json.load(open(RD + '/clinical.json'))['fi_tertile_cuts']
rows = ['SOFA ≥8', 'SOFA 4–7', 'SOFA 0–3']; cols = ['Low', 'Middle', 'High']
xl = [f'Low\n≤{c[0]:.2f}', f'Middle\n{c[0]:.2f}–{c[1]:.2f}', f'High\n>{c[1]:.2f}']
fig = plt.figure(figsize=(7.2, 6.4)); gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.72], left=0.11, right=0.90, bottom=0.08, top=0.93, wspace=0.12, hspace=0.62)
axs = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
for k, (ax, co, ttl) in enumerate(zip(axs, ['MIMIC-IV', 'eICU'], ['MIMIC-IV development', 'eICU external validation'])):
    d = g[g.cohort == co].set_index(['sofa', 'fi_tertile'])
    Z = np.array([[d.loc[(r, cc), 'mortality'] for cc in cols] for r in rows])
    im = ax.imshow(Z * 100, cmap='OrRd', vmin=0, vmax=45, aspect='auto')
    for i, r in enumerate(rows):
        for j, cc in enumerate(cols):
            q = d.loc[(r, cc)]; col = 'white' if q.mortality > 0.28 else 'black'
            ax.text(j, i - 0.12, f'{q.mortality * 100:.1f}%', ha='center', va='center', fontsize=11, fontweight='bold', color=col)
            ax.text(j, i + 0.22, f'{int(q.deaths)}/{int(q.n)}', ha='center', va='center', fontsize=8, color=col)
    ax.set_xticks(range(3)); ax.set_xticklabels(xl, fontsize=8.5); ax.set_yticks(range(3)); ax.set_yticklabels(rows if k == 0 else [])
    ax.set_xlabel('FI-Lab group (MIMIC-IV tertile cut-points)', fontsize=9); ax.set_title(ttl, fontweight='bold', fontsize=11, loc='left')
    ax.text(-0.14 if k == 0 else -0.04, 1.08, 'ab'[k], transform=ax.transAxes, fontsize=14, fontweight='bold')
    ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
p1 = axs[1].get_position(); cax = fig.add_axes([0.92, p1.y0, 0.015, p1.height]); cb = fig.colorbar(im, cax=cax); cb.set_label('Observed mortality, %')
o = pd.read_csv(RD + '/within_sofa_or.csv'); hj = json.load(open(RD + '/heterogeneity.json'))
ax = fig.add_subplot(gs[1, :]); CM, CE = '#2B6A94', '#B8642E'
pos = ax.get_position(); ax.set_position([pos.x0, pos.y0, 0.60, pos.height])
import matplotlib.transforms as mt
BT = mt.blended_transform_factory(fig.transFigure, ax.transData)
cats = ['SOFA 0–3', 'SOFA 4–7', 'SOFA ≥8']
for j, cat in enumerate(cats):
    for k, (co, col, dy) in enumerate([('MIMIC-IV', CM, 0.14), ('eICU', CE, -0.14)]):
        r = o[(o.cohort == co) & (o.sofa == cat)].iloc[0]; yy = len(cats) - 1 - j + dy
        ax.plot([r.lo, r.hi], [yy, yy], color=col, lw=1.6); ax.plot(r.OR, yy, 'o', color=col, ms=5, label=co if j == 0 else None)
        ax.text(0.74, yy, f'{co}: {r.OR:.2f} ({r.lo:.2f}–{r.hi:.2f})', fontsize=8.5, va='center', color=col, transform=BT, clip_on=False)
ax.axvline(1, color='#888', ls='--', lw=0.9); ax.set_xscale('log'); ax.set_xlim(0.9, 2.6)
ax.set_xticks([1, 1.25, 1.5, 2, 2.5]); ax.set_xticklabels(['1', '1.25', '1.5', '2', '2.5']); ax.minorticks_off()
ax.set_yticks(range(len(cats))); ax.set_yticklabels(cats[::-1]); ax.set_ylim(-0.6, len(cats) - 0.4)
ax.set_xlabel('Adjusted OR per 0.1 increase in FI-Lab (95% CI)')
for s_ in ['top', 'right']: ax.spines[s_].set_visible(False)
ax.legend(frameon=False, fontsize=8.5, loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=2)
pm, pe = hj['MIMIC-IV']['interaction_P'], hj['eICU']['interaction_P']
ax.set_title('Adjusted association of FI-Lab within SOFA categories', fontsize=10.5, fontweight='bold', loc='left', pad=22)
ax.text(-0.14, 1.12, 'c', transform=ax.transAxes, fontsize=14, fontweight='bold')
fig.savefig(FD + '/Figure4_risk_grid.png'); print('ok')
