# 脚本目录与执行顺序

以下按 pipeline 实际依赖顺序列出；6 与 6V 为互斥的插补分支。MR 为独立流程，不由 NHANES pipeline 自动启动。完整科学执行应另行安排新工作区。

|顺序|脚本（scripts/ 下）|输入|输出与用途|
|---|---|---|---|
|1|nhanes/01_build_cohort.R|public DEMO/BMX/CBC/SMQ/DPQ/MCQ/DIQ/BPQ/ALQ/RXQ_RX|work/cohort/C1_rebuilt_aligned_v3.rds; flow counts；sample 33174 and 2902 cases|
|2|nhanes/02_build_full_adult_design.R|full DEMO reconstruction|work/cohort/C2_full_adult_design.rds；full adult base 42143, degf109|
|3|nhanes/03_add_public_modules.R|cohort; public cotinine/CRP/SLQ/ALQ/drug dictionaries|work/cohort/C3_C1_new_variables_LOD_confirmed_pending_alcohol_antidepressants.rds；covariate inputs|
|4|nhanes/04_clinical_and_analyses.R|C2/C3|work/clinical/analytic_v2_antidepressant_pending.rds; clinical model CSVs；S1, CRP and age sensitivities|
|5|nhanes/15_run_D2_1.R|clinical covariates; RXQ/drug classification|work/analytic/analytic_v3.rds；antidepressant and complete analytic data|
|6|nhanes/14_run_sequence.R|analytic_v3; full adult design|work/imputation/imputations_standardized.rds; Rubin_MI_CC_comparison.csv；primary Model1–3 and Wald contrasts|
|6V|nhanes/verify_primary_from_frozen_imputations.R|rebuilt analytic_v3; historical frozen imputations|work/mi_validation/Rubin_MI_CC_comparison.csv；180 MI fits and9 CC fits;420 numerical comparisons|
|7|nhanes/06_prepare_supplement_inputs.R|C1 rebuilt|work/cohort/supplement_input.rds；logs/factor levels and PHQ15 field for supplementary models|
|8|nhanes/06_run_B1_2_to_6.R|supplement input;fulladult|work/symptoms_base/*.csv；S3/S4/S5 symptom and cycle analyses|
|9|nhanes/dependency_run_B2_1_and_audit.R<br>nhanes/07_run_B2_remaining.R|supplement input; B1;fulladult|work/symptoms_extended/*.csv；S3/S4 adjustments and symptom dimensions|
|10|nhanes/08_run_B3_and_ALQ_audit.R|supplement input;fulladult|work/thresholds/*.csv；PHQ thresholds|
|11|nhanes/09_run_B4.R<br>nhanes/dependency_run_B5.R<br>nhanes/dependency_run_B6.R|supplement input; public ALQ;fulladult|work/alcohol;work/education;work/inflammation；alcohol/socioeconomic/clinical sensitivities|
|12|nhanes/16_run_D2_3.R<br>nhanes/17_run_E1.R|analytic_v3; fulladult|work/analytic;work/descriptive；Table1; SII and sensitivities|
|13|nhanes/dependency_check_flow.R<br>nhanes/21_export_final.R<br>nhanes/22_export_remaining.R|prior results; fulladult|work/export/*.csv；Tables1/2; aggregate export;flow|
|14|nhanes/19_draw_figures.R<br>nhanes/20_draw_forest_publication.R|flow and primary CSVs|work/export/Figure1_flow.*; Figure2_forest_publication.*；Figures1/2|
|15|nhanes/25_run_G1.R|analytic_v3;fulladult|work/covariates/*.csv；S2 covariate attenuation; smoking; linear counts|
|16|nhanes/27_run_G2.R<br>nhanes/26_finalize_G2_v2.R|analytic_v3;fulladult;primary CSVs|work/discrimination/*.csv；AUC500 bootstrap;E-values|
|17|nhanes/28_run_G3.R<br>nhanes/29_run_G3_remaining.R<br>nhanes/30_run_G3b.R|analytic_v3;fulladult|work/shapley;work/heterogeneity；Shapley200 bootstrap;192 specs;indices;FiguresS1/S2|
|18|nhanes/32_run_H10_v3.R|analytic_v3;fulladult;G3 subset coefficients|work/attenuation/*.csv；attenuation1000 bootstrap|
|19|nhanes/36_run_phaseH7.R|analytic_v3;fulladult;imputations+scaling|work/race/*.csv；S5 race|
|20|nhanes/37_prepare_phaseH7b.R<br>nhanes/38_run_phaseH7b.R|public MCQ/KIQ_U/PAQ/SLQ;imputations;analytic_v3|work/sleep/*.csv；S4 additional covariates and sleep|
|21|nhanes/05_joint_splines.R<br>render_figure3.R|analytic_v3;fulladult;knots/grid/rug|work/spline/*.csv; Figure3.*；joint Figure3 and S8|
|22|nhanes/sex_interactions.R<br>nhanes/sex_total_score.R|analytic_v3;fulladult;imputations+scaling|work/sex;work/sex_total；revised S5 scalar and joint tests|
|23|export_revised_tables.py|output/work/sex; output/work/sex_total; output/work/spline (after their jobs)|workspace/output/tables_revised/S5_*.csv and S8_*.csv；revised S5/S8 direct rounding|
|MR1|mr/upstream_online_MR_UNVALIDATED.R|PGC file;OpenGWAS API;authentication;LD reference|work/mr_online；upstream instruments and MR|
|MR2|mr/refit_harmonized_MR.R|authorized harmonized UVMR/MVMR caches|work/mr_refit；MR model validation|
|MR3|mr/01_export_uvmr_diagnostics_final.R<br>mr/04_mvmr_final_export.R|historical authorized workspace/results/harmonized objects|output/work/mr_export only; inputs/gwas_cache is read-only；S6 diagnostic/strength lineage|
|MR4|mr/05_本地数据跑共定位分析.R|authorized NEU VCF+tbi;PGC file;harmonized N instruments|work/coloc/coloc_results_unique_regions_N563946.csv；S7 48regions|

## 通用工具与保留理由

- `init_workspace.py`：安全初始化、导入缓存与哈希记录。
- `run_job.py`：单个 R 入口、独占日志与工作区管理。
- `run_pipeline.py`：NHANES 依赖顺序；支持 --plan、--frozen-imputations；最终 S5/S8 导出和 Figure3 渲染均包含。
- `export_revised_tables.py`：从已拟合聚合 CSV 导出 S5/S8，直接舍入，无模型拟合。
- `render_figure3.R`：显式安全设备；默认可读取锁定聚合曲线，--curve-input 指向受控输出。
- `check_numeric_csv.py`：数值核对工具，使用前查看 --help。
- `lib/safety.R`、`lib/safety.py`：共有 I/O、配置、禁止覆盖及设备保护；新增归档父目录保护。
- `lib/domain_design.R`：全成人设计子总体实现。
- `environment/check_versions.R`：软件依赖实测版本检查。

`dependency_*.R` 名称不是废弃标记：它们是实际执行链的必要依赖，不可因名称删除。`26_finalize_G2_v2.R`、`32_run_H10_v3.R` 等虽含版本号，保留依据是已核验依赖和结果，不是文件名。

历史受限性别交互 CSV 仅作审计比较，无对应废弃分析入口重新加入本包。上游在线 MR 保留以完整记录方法，但明确 UNVALIDATED，不是正式数值重选入口。
