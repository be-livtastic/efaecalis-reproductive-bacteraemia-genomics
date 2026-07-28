FROM condaforge/miniforge3:24.11.3-0

LABEL org.opencontainers.image.title="E. faecalis ten-locus phylogeny"
LABEL org.opencontainers.image.description="Reproducible extended MLST-style whole-CDS analysis; not a core-genome pipeline"

WORKDIR /project
COPY environment/environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && mamba clean --all --yes

ENV PATH=/opt/conda/envs/efaecalis_phylogeny/bin:/opt/conda/bin:${PATH}
ENTRYPOINT ["/bin/bash", "scripts/06_phylogenomics/run_10_locus_phylogeny.sh"]
CMD ["--help"]
