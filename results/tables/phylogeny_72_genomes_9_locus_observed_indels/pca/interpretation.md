# Nine-locus alignment PCA

## Step 1: input validation

The PCA uses the completed 72-genome, 11,341-nt concatenated alignment. FASTA
identifiers are unique and match the 72 metadata records exactly. The source
alignment is read-only and is not modified.

## Step 2: feature construction

There are 320 variable alignment positions. Requiring each encoded alternative
allele to occur in at least two genomes produces 253 standardised
observed-allele features. Gap or ambiguous calls are treated as missing in the
numeric matrix and assigned the feature mean after centring (therefore zero);
this is a computational PCA convention and does not impute bases into any
sequence. The feature table retains locus, position, alleles, call counts and
missingness.

## Step 3: decomposition

SVD gives PC1 17.69%, PC2 13.36%, PC3 11.64%, PC4 10.58%, and PC5 9.51% of
the retained feature variance. PC1-PC2 together explain 31.05%; PCs 1-5 explain
62.77%. No single two-dimensional projection captures the complete structure.

yqiL contributes 45.8%, 34.5%, and 33.9% of the squared loadings on PCs 1-3,
respectively. PC2 also receives substantial contributions from gdh (20.1%) and
pyrC (19.1%). The main ordination structure is therefore strongly locus
weighted and should be interpreted as multilocus genetic structure, not as an
independent confirmation of the tree.

## Step 4: biological interpretation

Reproductive and bacteraemia isolates overlap in PCA space. Their PC1 means are
-1.75 and 0.42, while both categories span overlapping ranges; PC2 is strongly
affected by `GCA_029011255.1` and a bacteraemia isolate
(`GCA_050506875.1`). The ordination does not support a simple separation by
sample category. Labels are observational dataset categories and PCA alone
does not establish phenotype association, transmission, or causality.
