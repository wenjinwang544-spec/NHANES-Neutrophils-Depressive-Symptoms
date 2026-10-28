# System requirements and unresolved environment pinning

- Python3.11 or newer: pathlib.resolve/is_relative_to, AST/path unit tests and launchers. Tested through the bundled Python used in the audit; no application-specific runtime is required for ordinary launchers.
- R4.5.2: recorded study version; survey4.5, mice3.19.0 and the package inventory. jsonlite2.0.0 is required before managed study execution. R must run with --vanilla.
- macOS arm64/Quartz is the recorded platform. The unchanged Figure1/2 study renderers explicitly use Quartz. Linux device portability has not been tested. mclapply(mc.cores=4) requires fork support; the corresponding Windows execution is not claimed.
- shasum with -a256 is required by H10. tools::sha256sum from the recorded R runtime is required by input manifests. File readers require haven/data.table/seqminer as listed.
- OpenGWAS upstream extraction requires authorized access and suitable credentials. The retained historical upstream script references IEUGWAS_TOKEN; verify the credential mechanism against the installed ieugwasr version before any authorized online run. The environment variable name alone does not prove authenticated access. No token was read or requested during this repair.
- ieugwasr1.1.0.9000 is a development build whose exact historical-analysis source commit is not established (the current installed RemoteSha is now listed in environment/INSTALL.md). package_versions.csv is an observed inventory, not a complete cross-platform installation lockfile. This gap remains explicit.

Safety tests use synthetic files and no clinical statistical fits. Formatting checks require python-docx only in the private review environment, not to run the scientific model pipeline. No current installed package or environment file was modified.


## Final archive update (2026-10-28)

See docs/VALIDATION_STATUS.md for Stage8 limited execution, README.md for current initialization, and docs/RELEASE_NOTES.md for this assembly. The archive adds automatic source-parent protection and decodes R --file spaces in managed headers; no scientific body or locked result is changed. Full pipeline/online MR remain outside the new execution claim. The scoped original-code MIT text is prepared; actual publication awaits final approval.
