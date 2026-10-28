# Final code ownership and license-scope review

## Prepared notice and scope

Copyright (c) 2026 Wenjin Wang. The complete, unmodified MIT permission/warranty text is in LICENSE. The name is the author's current specified spelling; final approval to publish is still withheld.

The MIT license in LICENSE applies only to original software code that Wenjin Wang is entitled to license: scripts/**, tests/**, environment/check_versions.R, and the original configuration config/storage_policy.json. It does not relicense third-party packages or their implementations, external data, figures, aggregate/reference CSVs, or other research assets. Those materials retain their respective rights and terms. No permission to redistribute GWAS inputs, participant-level data or private caches is granted. The package has not been published; this prepared candidate awaits the author's final release approval.

## Examined source and findings

All52 R/Python source files were reviewed, including 49 scripts/libraries, two tests and environment/check_versions.R; the configuration was checked separately. The 41 scientific entrypoints and all other executable sources are byte-identical to the accepted147-file package. Imports/calls to survey, mice, Hmisc, MVMR, TwoSampleMR and other packages are not bundled copies of their implementations. Third-party packages retain their listed licenses. No named other substantive human contributor, embedded external license header or identifiable substantial copied implementation was found.

The one positive provenance question was scripts/mr/upstream_online_MR_UNVALIDATED.R: the selection index described extracting a block from a monolithic script that also contained pasted analyses. The author directly confirmed that this selected MR block was original and did not copy substantive third-party code. That question is closed on the author's declaration; it is not inferred from the absence of a copyright header.

Source extraction/adaptation from the author's own historical research files is documented by provenance/script_selection.csv. Safety wrappers/launchers/tests were produced during the AI-assisted audit. No separate named contributor or external copied implementation was identified for those additions; AI assistance does not by itself establish every possible copyright interest. The grant is limited to rights the author is entitled to license, and does not erase third-party rights if a later claim or new source is found.

## E-value reference: method citation versus code copying

scripts/nhanes/26_finalize_G2_v2.R:17 contains an implementation-reference URL to EValue and a paper citation. scripts/nhanes/27_run_G2.R:60–63 implements the short published risk-ratio E-value formula and records that methodological source. The reference EValue package is GPL-2, not MIT. No substantial EValue function body, documentation example, package source or dataset is included; use of the same mathematical formula is not evidence that the full GPL implementation has been copied. Preserve the existing method citation. If the author later identifies actual copied EValue implementation, it must be reviewed under its GPL terms rather than assumed MIT.

Sources checked: https://louisahsmith.github.io/evalue/reference/evalues.OR.html ; https://louisahsmith.github.io/evalue/ ; https://raw.githubusercontent.com/cran/EValue/master/DESCRIPTION ; https://raw.githubusercontent.com/cran/EValue/master/R/EValue.R .

## Figures, data and other rights

Figure3 PDF, all reference/aggregate CSVs and external data are outside the MIT scope. The PNG with an embedded Apple ICC profile remains excluded. No attribution was stripped, no image was edited and no model or figure was regenerated. External GWAS use terms and dependency licenses are separate from the original-code grant.

No identifiable unresolved third-party source/authorization item remains after the author's MR confirmation. This does not certify that every unmarked historical snippet is original; future source discoveries require a new review. The grant covers only the author's rights and does not promise rights in other people's work.

MIT reference: https://opensource.org/license/mit . No third-party GPL source is being relicensed as MIT. This review does not determine every legal issue concerning a combined distribution of external dependencies.
