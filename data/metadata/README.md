# Sample metadata

`sample_metadata_source_72.tsv` is the canonical public sample catalog. It was
exported from the selected `reproductive samples` and `Bacteraemia samples`
sheets in the local source workbook, with whitespace normalized and the country
typo `United Kingdon` corrected to `United Kingdom`. It contains no personal
paths or workbook row numbers. Collection years available in the workbook's
public metadata sheet or retained NCBI assembly reports were added for 15
samples; unavailable years remain blank.

`create_supplementary_sample_table.R` validates that the catalog contains the
same 72 assembly/category pairs as the accession manifest, requires unique and
well-formed assembly, BioSample and strain identifiers, and enforces the 14/58
group counts. It then regenerates the analysis metadata, identifier-QC report
and manuscript Supplementary Table S1.

The workbook-derived assembly–BioSample–strain mappings exactly matched the
older archived metadata export. The 15 locally retained NCBI assembly reports
also matched for every identifier they contained; the other mappings retain
their workbook/public-metadata provenance.
