# FI-Lab AMI external validation — analysis code (manuscript v0.60)

Code for *External validation of a laboratory-based frailty index beyond SOFA for 30-day in-hospital mortality in critically ill patients with acute myocardial infarction*.

A harmonized 27-item laboratory frailty index (FI-Lab) is scored from first measurements in MIMIC-IV and eICU. Logistic models with clinical covariates and SOFA, with and without FI-Lab, are developed in MIMIC-IV patients alive and hospitalized beyond 24 hours after ICU admission. The models are then applied unchanged to the corresponding eICU patients.

Reference ranges in `spec/item_crosswalk_and_thresholds.csv` follow the published 33-item FI-Lab of Ding et al. (BMC Infect Dis 2024;24:573; doi:10.1186/s12879-024-09430-w); the male creatinine range (0.74–1.35 mg/dL) is the Mayo Clinic Laboratories adult range, because the source table lists 3.5–5.2 mg/dL.

The repository contains code, specifications and aggregate reference results only. It contains no patient-level data.

## Requirements

- Python ≥3.10 with `pip install -r requirements.txt` (numpy, pandas, openpyxl, matplotlib).
- R ≥4.1 with `pROC`, `ggplot2` and `patchwork`, for the main analysis and the component-threshold analysis (stage 2).
- Credentialed PhysioNet access to MIMIC-IV v3.1 and the eICU Collaborative Research Database. Inputs are three exports produced by the authors' SQL extraction:

| Argument | Export | Content |
|---|---|---|
| `--mimic` | `mimic.xlsx` | MIMIC-IV AMI ICU records: first and window-extreme laboratory and vital-sign values, SOFA, APS III, comorbidities, outcome |
| `--eicu` | `eicu_max及min.xlsx` | eICU AMI ICU records: window-extreme values, SOFA, APACHE IV, comorbidities, outcome |
| `--eicu-first` | `eicu首测值xlsx.xlsx` | eICU first measurements, in the same row order as `--eicu` |

## Run

```sh
python run.py --mimic /private/mimic.xlsx --eicu /private/eicu_max及min.xlsx \
  --eicu-first /private/eicu首测值xlsx.xlsx --output /private/filab_run
```

The output directory must be new and must be outside the repository, because it contains patient-level data. Options:

- `--skip-r` runs only the Python stages. In this mode, the figures use `reference_results/main` for the R-derived inputs.
- `--skip-figures` skips all plotting.
- `--rscript PATH` selects the R executable.

The pipeline runs in five stages:

1. **`scripts/prepare.py`** prepares the cohorts:
   - Scores FI-Lab. Abnormal components are divided by observed components, and at least 18 of 27 must be observed.
   - Selects one admission per patient. In MIMIC-IV this is the earliest admission. In eICU it is a random admission, using an R-compatible Mersenne-Twister with seed 20261022.
   - Derives the outcome and covariates, and applies the 24-hour survivor rule.
   - Writes `data/{M,E}_design.csv`, `data/cohorts.pkl` and `results/cohort_flow.json`.
2. **`scripts/main/*.R`** runs the main analysis: frozen external validation with 1,000 paired bootstrap samples, internal optimism correction (200 refits), the DeLong test, baseline summaries, split-half recalibration and decision curves. **`scripts/additional/threshold.R`** (after `threshold_export.py`) repeats the frozen-model comparison at minimum observed-component thresholds of 12–18 and compares patients by component coverage (Tables S7–S8).
3. **`scripts/additional/analyse.py`** and **`renal_sofa.py`** run the additional analyses: creatinine-based renal SOFA (Table 2), associations and splines, external sensitivity analyses, the APACHE IV and APS III comparisons, overlap with severity scores, and component summaries.
4. **`scripts/additional/figures.py`** and **`figure5.py`** draw Figures 1–4 and Supplementary Figures S1–S4. **`clinical.py`** produces the risk grid, 10% threshold classification, standardized risks and the three-variable model (Figure 4a–b, Table 3); **`heterogeneity.py`** produces FI-Lab associations within SOFA categories with interaction tests and region-stratified validation in eICU (Figure 4c, region rows of Figure 3, Tables S9 and S13).
5. **`tests/compare_reference.py`** compares all aggregate results with `reference_results/` (tolerance 1e-6).

## Manuscript map

| Manuscript item | Produced by | Reference file |
|---|---|---|
| Figure 1 (study design and cohort flow) | prepare.py, figures.py | `reference_results/cohort_flow.json` |
| Table 1 | describe.R, analyse.py (`overlap`) | `main/baseline.csv`, `additional/overlap.json` |
| Table 2 (exported SOFA), Figure 2 | analysis.R, analyse.py (`eicu_pred_main.csv`), figures.py | `main/external_validation.csv`, `main/internal_validation.csv`, `main/external_auc_tests.csv`, `additional/renal_sofa.csv` |
| Table 2 (creatinine-based renal SOFA) | renal_sofa.py | `additional/renal_sofa.csv` |
| Table S5, Supplementary Figure S3 (splines) | analyse.py, figures.py | `additional/associations.json` |
| Figure 3 | analyse.py, renal_sofa.py, heterogeneity.py, figures.py | `additional/external_sensitivity.csv`, `additional/renal_sofa.csv`, `additional/within_severity.json`, `additional/eicu_region_validation.csv` |
| Figure 4a–b, Table 3 | clinical.py, figure5.py | `additional/risk_grid.csv`, `additional/clinical.json` |
| Figure 4c, Figure 3 (region rows), Tables S9 and S13 | heterogeneity.py | `additional/within_sofa_or.csv`, `additional/eicu_region_validation.csv`, `additional/heterogeneity.json` |
| Table S3 (coefficients) | analysis.R | `frozen_coefficients.csv` (run output) |
| Supplementary Figure S1 (FI-Lab distribution) | analyse.py, figures.py | `fi_values.csv` (run output) |
| Table S4, Figure S2 | analyse.py | `additional/components.csv` |
| Table S6 | analyse.py | `additional/external_sensitivity.csv` |
| Tables S7–S8 (component thresholds) | threshold_export.py, threshold.R | `additional/threshold_performance.csv`, `additional/included_excluded.csv` |
| Table S10 | analyse.py | `additional/within_severity.json`, `additional/overlap.json` |
| Tables S11–S12, Figure S4 | recalibration_dca.R | `main/recalibration_survivors.csv`, `main/decision_curves_survivors.csv` |

The exported-SOFA estimates in the manuscript (Table 2 upper block, abstract) come from R; the creatinine-based block of Table 2 comes from renal_sofa.py. The additional analyses use NumPy and reproduce the same cohorts, coefficients and point estimates. Their bootstrap intervals for the main comparison differ slightly from the R intervals, because the two use different random-number streams. Table S7 reports the ≥18 rows from Table 2.

## Data protection

Do not commit or share the contents of the output directory. This includes the design matrices, predictions, bootstrap draws and `cohorts.pkl`. MIMIC-IV and eICU data may only be used under their PhysioNet data use agreements.
