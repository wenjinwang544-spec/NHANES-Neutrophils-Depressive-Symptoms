# Storage configuration

storage_policy.json is source policy: protect_source_parent=true automatically protects this repository's containing directory. It is not a workspace config or security marker. Actual local absolute paths are supplied to scripts/init_workspace.py; its generated private run_config.json lives only in the external workspace. Do not put account paths, credentials or imported data here.
