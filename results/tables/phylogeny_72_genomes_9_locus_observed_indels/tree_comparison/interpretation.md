# Comparison of the three completed phylogenies

The comparison uses SH-aLRT >=80% and ultrafast bootstrap >=95% as the joint
support threshold. A major split has at least four taxa on its smaller side.
Trees are compared on their 71 shared genomes so removal of
`GCA_029011395.1` is not mistaken for topological disagreement.

## Supported major structure

Each tree contains six supported major splits. Five are identical across all
three trees. One additional 23-genome split is shared by the primary and
71-genome nine-locus trees but not the seven-locus tree; the seven-locus tree
instead has one unique supported four-genome split. Thus the broad supported
structure is largely reproducible, while a limited part of the internal
resolution depends on retaining yqiL and recA.

## Reproductive-isolate placement

The reproductive isolates are dispersed rather than forming a single supported
source-specific clade. In the primary tree, 12 of 14 have a bacteraemia isolate
as their nearest neighbour. The exception is the strongly supported, nearly
identical pair `GCA_029011395.1`/`GCA_029011745.1`. After the former is removed,
the latter's nearest neighbour necessarily changes to a bacteraemia isolate.

Most placements are stable. The seven-locus tree changes the nearest neighbour
for `GCF_049262565.1`, `GCF_053117665.1`, `GCF_053117695.1`, and
`GCF_906464925.1`; two of those changes create a reproductive-reproductive
nearest-neighbour pair. These are local rearrangements and should not be
described as source-group separation.

## Branch-length behaviour

The primary and 71-genome nine-locus trees have nearly identical total lengths
(0.087507 and 0.087586 substitutions/site). The seven-locus tree is shorter
(0.076416), as expected when two informative loci are removed, so raw lengths
should not be compared as if based on the same character set.

`GCA_050472755.1` has the longest terminal branch in every tree. Other repeated
long-terminal taxa include `GCF_053117695.1`, `GCA_050506775.1`, and
`GCA_029011255.1`. Removing `GCA_029011395.1` makes its close partner
`GCA_029011745.1` terminally long; this is a pruning effect, not evidence that
the retained sequence evolved suddenly. Many terminal branches are zero or
near zero because several genomes have identical multilocus sequences, making
median/MAD outlier rules unstable. Rankings and absolute branch lengths are
therefore reported instead.
