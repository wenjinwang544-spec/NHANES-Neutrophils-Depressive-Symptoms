# Required input contracts

## NHANES public modules

Provider: https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/ . Seven cycles2005–2006 through2017–2018, suffixesD–J. Use `--raw-nhanes` to configure the existing public-data root. Readers use cycle/component folders and case-insensitive filenames:

```
2005-2006/Demographics/demo_d.xpt
2005-2006/Examination/bmx_d.xpt
2005-2006/Laboratory/cbc_d.xpt
2005-2006/Questionnaire/dpq_d.xpt
```

Required core modules: DEMO, BMX, CBC, SMQ, DPQ, MCQ, DIQ, BPQ, ALQ, RXQ_RX. Additional inputs: COT in2005–2006 and2013–2018, COTNAL in2007–2012; CRP in2005–2010 and HSCRP in2015–2018, with unavailable cycles retained as missing according to the unchanged scientific code. SLQ and RXQ_DRUG are read across cycles; RXQ_DRUG is the unsuffixed dictionary file. Sleep/chronic-condition preparation additionally reads MCQ, KIQ_U, PAQ and SLQ for all seven cycles. Its new manifest must contain28 actual cycle-module files, nonempty SHA256 and sizes, and must fail if empty or incomplete.

The initializer does not download or alter public data. No participant-level NHANES data are shipped in this repository. Scientific counts and treatment of missingness are unchanged.

## Frozen imputation validation

Import `imputations_standardized.rds` with init_workspace.py. It is the historical mice object, not a newly reconstructed dataset. The frozen validation script derives/writes its scaling and model outputs into the output tree; it never alters the imported object. For individual downstream interaction jobs, `imputation_scaling.csv` may also be imported if not already generated in the output tree. The object schema/row alignment is verified by the unchanged scientific code; original scientific generation and seed are documented separately.

## Genetic inputs, never bundled

Use repeated `--import-cache filename=/authorized/path` during initialization. Imported analysis inputs live in `inputs/gwas_cache`, separate from `output`.

- `refit_harmonized_MR.R`: harmonised_list.rds keyed by the three exposure labels with TwoSampleMR harmonized columns; mvmr_harmonised_final.rds with the retained multivariable harmonized object. These are frozen-object validation inputs, not independent instrument extraction.
- `04_mvmr_final_export.R`: .RData containing mvmr_harmonised and mvmr_results.csv with id.exposure, b, se, pval and nsnp. Its exported results and retained harmonized output go only to output/work/mr_export; the input cache stays unchanged.
- `01_export_uvmr_diagnostics_final.R`: .RData containing harmonised_list, heterogeneity_list, pleiotropy_list, steiger_list and uvmr_results_list; harmonised_list.rds; locked uvmr_results_full.csv; and the historical mr_noUKBB修正版.R text needed for provenance checks. No historical workspace is distributed here.
- Local colocalization: ieu-b-34.vcf.gz and matching .tbi, pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11_formatted.tsv.gz, and harmonised_list.rds. `data/gwas/` logical reads resolve to the explicit imported cache, not an input/output symlink. Preserve exact provider versions and checksums.
- Upstream online MR: local PGC outcome file plus authorized OpenGWAS exposure/LD services and credentials; retained code remains UNVALIDATED. Re-downloads may not reproduce historical instruments. Exact historical random seeds and the development ieugwasr commit remain unresolved.

All inputs remain subject to provider terms. See DATA_SOURCES_AND_RIGHTS.md; no code license grants permission to redistribute these files. Do not place data, .RData, local run configs, tokens or manuscripts in a Git commit.
