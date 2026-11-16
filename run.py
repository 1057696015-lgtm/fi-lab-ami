"""Reproduce all analyses for the FI-Lab AMI external validation manuscript (v0.60).

python run.py --mimic mimic.xlsx --eicu eicu_max及min.xlsx --eicu-first eicu首测值xlsx.xlsx --output /private/run

Stages: (1) prepare cohorts [Python]; (2) main analysis and component-threshold analysis [R]; (3) additional analyses and figures [Python].
--skip-r runs stages 1 and 3 only, using reference_results/main for the R-derived inputs of the figures.
The output directory must be new and outside this repository: it contains patient-level data.
"""
import argparse, shutil, subprocess, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for n in ['mimic', 'eicu', 'eicu-first', 'output']: p.add_argument('--' + n, required=True)
p.add_argument('--rscript', default='Rscript'); p.add_argument('--skip-r', action='store_true'); p.add_argument('--skip-figures', action='store_true')
a = p.parse_args(); out = Path(a.output).expanduser().resolve()
if out == ROOT or ROOT in out.parents: p.error('Use an output directory outside the repository (outputs are patient-level).')
if out.exists(): p.error('Output directory already exists; use a new one.')
for n in ['mimic', 'eicu', 'eicu_first']:
    if not Path(getattr(a, n)).expanduser().is_file(): p.error('Missing input: ' + n)
for sub in ['data', 'results/main', 'results/additional', 'figures', 'scripts', 'logs']: (out / sub).mkdir(parents=True, exist_ok=True)

def run(name, cmd):
    print('>>', name, flush=True)
    with (out / 'logs' / f'{name}.log').open('w') as log:
        r = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if r.returncode: print((out / 'logs' / f'{name}.log').read_text()[-3000:]); sys.exit(f'{name} failed')

py = [sys.executable, '-I']
run('01_prepare', py + [str(ROOT / 'scripts/prepare.py'), '--mimic', a.mimic, '--eicu', a.eicu, '--eicu-first', a.eicu_first, '--spec', str(ROOT / 'spec'), '--output', str(out)])

main_results = ROOT / 'reference_results/main'
if not a.skip_r:
    if not shutil.which(a.rscript): p.error('Rscript not found (install R with pROC, ggplot2, patchwork) or use --skip-r')
    for f in (ROOT / 'scripts/main').glob('*.R'): shutil.copy2(f, out / 'scripts' / f.name)
    (out / 'results').mkdir(exist_ok=True)
    for name in ['analysis.R', 'describe.R', 'recalibration_dca.R']:
        run('02_' + name, [a.rscript, str(out / 'scripts' / name)])
    for f in (out / 'results').glob('*.*'): shutil.move(str(f), out / 'results/main' / f.name)
    main_results = out / 'results/main'
    run('02_threshold_export', py + [str(ROOT / 'scripts/additional/threshold_export.py'), str(out / 'data/cohorts.pkl'), str(out / 'data')])
    run('02_threshold', [a.rscript, str(ROOT / 'scripts/additional/threshold.R'), str(out / 'data'), str(out / 'results/additional')])

run('03_additional', py + [str(ROOT / 'scripts/additional/analyse.py'), str(out / 'data/cohorts.pkl'), str(out / 'results/additional')])
run('03b_renal_sofa', py + [str(ROOT / 'scripts/additional/renal_sofa.py'), str(out / 'data/cohorts.pkl'), str(out / 'results/additional')])
run('03c_clinical', py + [str(ROOT / 'scripts/additional/clinical.py'), str(out / 'data/cohorts.pkl'), str(out / 'results/additional')])
run('03d_heterogeneity', py + [str(ROOT / 'scripts/additional/heterogeneity.py'), str(out / 'data/cohorts.pkl'), str(out / 'results/additional')])
if not a.skip_figures:
    run('04b_figure5', py + [str(ROOT / 'scripts/additional/figure5.py'), str(out / 'results/additional'), str(out / 'figures')])
    run('04_figures', py + [str(ROOT / 'scripts/additional/figures.py'), str(out / 'results/additional'), str(main_results),
                             str(out / 'data/cohorts.pkl'), str(out / 'figures'), str(out / 'results/cohort_flow.json')])
run('05_compare_reference', py + [str(ROOT / 'tests/compare_reference.py'), str(out)])
print('Completed:', out)
