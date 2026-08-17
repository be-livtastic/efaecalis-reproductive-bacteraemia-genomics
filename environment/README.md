# Software environment and version records

`environment.yml` is the canonical reproducible environment specification. Every direct Conda and pip dependency is pinned to an explicit version. The file was validated with a non-mutating Mamba dry run against the installed `efaecalis_phylogeny` environment.

`software_versions.tsv` mirrors every direct package and version in `environment.yml`. Update these two files together whenever a dependency pin changes.

`recorded_analysis_versions.tsv`, `amrfinder_r_session_info.txt`, and `amrfinderplus_version.tsv` are historical analysis provenance. Their versions can legitimately differ from the current environment specification and must not be rewritten to imply that an already completed analysis used a newer package. The historical AMR R session was run under Windows with R 4.5.2, whereas the current pinned Linux environment uses R 4.5.3.

Samtools 0.1.19 is intentionally pinned because Prokka 1.15.6 currently resolves through the legacy `perl-bio-samtools` dependency. A future upgrade to modern Samtools requires separating Prokka into another environment or container and validating both workflows before changing this pin.

AMRFinderPlus database releases are data resources rather than Conda package versions. The database version used for the completed analysis remains recorded in `amrfinderplus_version.tsv` and `recorded_analysis_versions.tsv`.
