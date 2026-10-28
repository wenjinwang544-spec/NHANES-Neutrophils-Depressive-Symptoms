# Environment setup

Use a separate software environment and library, outside this source archive. Python >=3.11 is sufficient for launchers, exporters and safety tests (standard library only). The audit used Python3.12.14. R4.5.2/macOS arm64 is the recorded study platform.

Inspect environment/package_versions.csv and environment/runtime.txt first. Most ordinary packages can be installed from CRAN into your separate library, for example:

```r
install.packages(c("survey", "mice", "haven", "dplyr", "tidyr", "data.table", "ggplot2", "jsonlite", "seqminer", "coloc"), repos="https://cloud.r-project.org")
```

This example installs available versions, not necessarily the recorded versions. Review and install matching archived/source versions where required, then compare packageVersion() with package_versions.csv. Do not silently accept changed statistical outputs after a package update. Read the remaining imported package names from docs/dependency_imports.csv.

TwoSampleMR0.6.29 and MVMR0.4.2 were retained locally; obtain the matching official upstream release source (MRCIEU/TwoSampleMR and WSpiller/MVMR), verify version/source and satisfy their dependencies before use. This archive does not include third-party package source. ieugwasr1.1.0.9000 is a development version whose exact historical commit was not recorded: a complete exact installation recipe is unavailable, and online upstream reconstruction remains UNVALIDATED. No renv lockfile is invented.

No installer runs automatically. Rscript must use --vanilla. `jsonlite` is required for path configuration before study code. `shasum -a256` is required by H10. Bootstrap fork workers and Quartz plots have only the recorded macOS validation. The package inventory is observed rather than an exhaustive clean-machine lockfile.

Do not place a private .Renviron in the repository. Historical online code references IEUGWAS_TOKEN; authentication details depend on the installed ieugwasr version and must be checked before any separately authorized online run. No token is required for offline frozen exports, and no token is embedded here.

## Current installation provenance observed during publication review

All153 entries in the package inventory match the current installed versions. R4.5.2/macOS arm64 and Python3.12.14 were rechecked without loading study models. Current installed DESCRIPTION files report these RemoteSha identifiers:

- ieugwasr1.1.0.9000: 8d347809d3eb9f00ed5ebc7861008401a000c900
- TwoSampleMR0.6.29: bd4d9d59431fc513722a2141d00552fa2d995a6f
- MVMR0.4.2: c7225eb61fd4843127cd9e005be6114dc4ae364d

These identify the installation observed now, not independently proven historical analysis commits. A clean-machine install at these commits and matching numerical reproduction have not been performed. The observed commit information narrows the installation gap but does not resolve the historical provenance gap or validate online MR. No packages were installed or updated.
