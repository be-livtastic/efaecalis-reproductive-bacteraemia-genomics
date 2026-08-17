# HLGR-associated genotype proxy decision

## Decision

For this dissertation, an accepted genomic detection of `aac(6')-Ie/aph(2'')-Ia` is labelled **HLGR-associated genotype proxy detected**. Absence from the validated functional AMR matrix is labelled **HLGR-associated genotype proxy not detected**.

These labels are not `HLGR+` and `HLGR-`. Phenotypic high-level gentamicin resistance requires antimicrobial susceptibility testing and cannot be confirmed from genomic detection alone.

## Analytical source

Proxy calls are read from the QC-passed AMR outputs. The MLST integration does not reinterpret raw AMRFinderPlus records. This preserves the upstream accepted-hit policy and its audit trail.

## Evidence and limitation

The bifunctional determinant is strongly associated with high-level gentamicin resistance in enterococci, including *E. faecalis*. However, gene-positive isolates without the corresponding high-level phenotype have been reported when the locus is disrupted. The dissertation therefore treats detection as a genotype proxy and preserves the AMRFinderPlus method, coverage, identity and contig supporting each call.

Primary supporting sources include:

- Daikos et al. (2003), DOI: <https://doi.org/10.1128/AAC.47.12.3950-3953.2003>
- Hsieh et al. (2021), DOI: <https://doi.org/10.1093/jac/dkab092>
- Chow et al. (2001), multiplex detection of enterococcal gentamicin-resistance genes: <https://pmc.ncbi.nlm.nih.gov/articles/PMC152526/>

## Future determinant review

The version-controlled configuration seeds a review watchlist with `aph(2'')-Ib`, `aph(2'')-Ic`, and `aph(2'')-Id`. A future accepted non-primary aminoglycoside determinant that is watchlisted or annotated as gentamicin-associated stops final proxy classification as `HLGR_proxy_review_required`. Its biological role must be reviewed before the proxy definition is expanded.
