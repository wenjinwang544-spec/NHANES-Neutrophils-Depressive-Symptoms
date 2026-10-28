# Run guide and input/output contracts

Start at repository root, configure a nonexistent absolute workspace outside the containing research directory and all inputs, and pass all historical roots explicitly. config/storage_policy.json automatically protects the source parent in Python and direct R entrypoints. Copying this archive elsewhere also protects its new containing parent; choose a workspace outside that parent. Do not remove this policy to fit an invalid path.

## Frozen NHANES branch

Initialize with --raw-nhanes plus `--import-cache imputations_standardized.rds=/authorized/object.rds`. Then --plan, check-only and, for a separately arranged analysis run, run_pipeline.py --frozen-imputations. The pipeline builds analytic/fulladult prerequisites before downstream jobs. Do not supply a different mice object merely to make a script run.

For an isolated downstream sleep/sex task, prerequisite work/analytic/analytic_v3.rds and work/cohort/C2_full_adult_design.rds must already be present in that managed workspace. Imported imputations_standardized.rds/imputation_scaling.csv have explicit read-only fallback mappings. Other logical work/* objects are not automatically remapped from the input cache. Stage8 staged authorized frozen copies through the exclusive output() API and verified hashes; no staging shortcut or symlink is embedded here. Missing prerequisites must stop. Never copy the Stage8 private workspace wholesale.

## MVMR-only frozen export

Use a new workspace (separate from any earlier MR export run):

```sh
python3 scripts/init_workspace.py --workspace /separate/new_mvmr_session \
  --protected-root /path/to/research_root \
  --protected-root /path/to/stage7 \
  --import-cache .RData=/authorized/final_workspace.RData \
  --import-cache mvmr_results.csv=/authorized/mvmr_results.csv
python3 scripts/run_job.py scripts/mr/04_mvmr_final_export.R --workspace /separate/new_mvmr_session --rscript /path/to/Rscript
```

.RData must contain mvmr_harmonised; CSV requires id.exposure,b,se,pval,nsnp. Script uses mv_multiple only to confirm the frozen object agrees with saved estimates, then exports strength and results to output/work/mr_export. Stage8 b/SE/P match is not an online instruments reconstruction.

For harmonized UVMR refit, separate new workspace imports harmonised_list.rds and mvmr_harmonised_final.rds. Historical stochastic SE/P are not authoritative replacements. UVMR diagnostic export needs additional historical workspace objects and text provenance files, listed in data/README.md; it was not rerun in Stage8.

Local colocalization imports the precise VCF, its .tbi index, named PGC file and harmonised_list.rds. Its Stage5 component result matched 48 regions byte-for-byte. External licensing and offline availability are independent requirements. No external genetic input is distributed.

## S5/S8 and figures

The full NHANES pipeline already schedules spline, sex_interactions.R and sex_total_score.R before export_revised_tables.py. Standalone export needs six aggregate CSVs in managed output/work/{sex,sex_total,spline}; the canonical reference copies live in supplementary/validated. These contain aggregate coefficients/curves only, not mice data. The exporter formats current full three-count tests, not the audit-only restricted sex CSV. It fails if inputs are missing or tables_revised already exists.

Standalone Figure3 uses the locked aggregate files unless --curve-input specifies a different validated source. It is not an automatic refit or new statistical result. Run through run_job.py so that outputs are explicit and exclusive. Figure1/2 are downstream of model/export prerequisites.

## Safe no-model checks

`tests/test_archive_layout.py` creates only synthetic fixtures under docs/validation, tests path guards and direct uninitialized entrypoints, preserves sanitized summary/failure evidence, and removes only its newly created synthetic fixture directory on success. It never loads study inputs or fits models. Existing real scientific workspaces remain untouched.

The inherited tests/test_output_safety.py is the Stage7 portable synthetic suite; it creates temporary fixtures outside source and may be used when such a temporary location has been authorized. It was not executed during archive assembly because this session permits writes only within this new folder.

All failure logs are evidence. Do not retry a real failed job into existing outputs: investigate, report, then obtain an appropriate new-run decision and a fresh workspace. Check-only verifies path setup, not scientific input completeness or model reproducibility.
