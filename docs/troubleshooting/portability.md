# Portability notes

- Run Bash bioinformatics tools under Linux or WSL unless their documentation
  explicitly supports the host platform.
- Activate the pinned environment before running a stage.
- Paths are resolved from each script or from `EFAECALIS_PROJECT_ROOT`; do not
  insert a personal absolute path into a committed script.
- Linux/WSL path and filename case matters. Use `AMRFinder`, not a mixed-case
  spelling that only works on Windows.
- Raw and interim outputs are local-only. Existing destinations trigger a safe
  failure; review them rather than deleting or overwriting them automatically.
- Use `config/paths.yml` locally only when overriding the example configuration.

