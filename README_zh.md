# FI-Lab AMI 外部验证分析代码（稿件 v0.60）

本包只含代码、规格文件和汇总参考结果，不含任何患者级数据。完整说明见 README.md。

FI-Lab 参考范围（spec/item_crosswalk_and_thresholds.csv）取自 Ding 等 2024 年 BMC Infectious Diseases 发表的 33 项 FI-Lab（doi:10.1186/s12879-024-09430-w）；原表男性肌酐范围为 3.5–5.2 mg/dL，不是成人参考范围，故改用 Mayo Clinic Laboratories 成人男性范围 0.74–1.35 mg/dL。

## 运行

    python run.py --mimic mimic.xlsx --eicu eicu_max及min.xlsx --eicu-first eicu首测值xlsx.xlsx --output /私有目录/filab_run

输出目录必须是新目录，且位于本仓库之外，因为里面有患者级中间数据。没有 R 环境时加 `--skip-r`，只运行 Python 部分；作图时 R 部分的输入改用 reference_results/main 中的参考结果。

## 流程

1. prepare.py：FI-Lab 评分（27 项中至少观测到 18 项）；每位患者选一次住院（MIMIC-IV 取最早一次，eICU 按种子 20261022 随机选取，与 R 结果逐行一致）；生成结局和协变量；按 24 小时存活规则筛选。
2. scripts/main/*.R：主分析，包括冻结模型外部验证、内部校正、DeLong 检验、再校准和决策曲线；scripts/additional/threshold.R（先运行 threshold_export.py）：最少观测项目数 12–18 的门槛分析及按覆盖程度的患者比较（Table S7、S8）。
3. analyse.py 和 renal_sofa.py：附加分析，包括按肌酐重算肾分项的 SOFA（Table 2 下半部分）、关联与样条、敏感性分析、APACHE IV/APS III 比较、各检验项目汇总。
4. figures.py、figure5.py：绘制 Figure 1–4（Figure 1 为研究设计与纳入流程，Figure 2 为外部验证 ROC、校准及内外部增量）和 Figure S1–S4（S1 FI-Lab 分布，S2 各项目观测，S3 样条曲线，S4 决策曲线）；clinical.py：风险网格、10% 阈值分类、标准化风险和三变量模型（Figure 4a–b、Table 3）；heterogeneity.py：SOFA 分层内 FI-Lab 关联及交互检验、eICU 分地区外部验证（Figure 4c、Figure 3 地区行、Table S9 与 S13）。
5. tests/compare_reference.py：将运行结果与参考结果核对，容差 1e-6。
