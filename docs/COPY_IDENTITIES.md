# Specific copy identity and protein-family similarity

Reports preserve a specific locus name separately from protein-database naming,
the user-selected comparison group and a locally editable copy label. A broad
NANOG hit must not rename a known NANOG-like locus to NANOG. Numbered paralogs
and `-like` suffixes are never removed to make names agree.

`genecontext/gene_identity.py` provides derived identity descriptions without
modifying scientific evidence or decisions. `gene_identities.tsv` exposes the
original identifiers, label, status, basis, protein-database names, discrepancies
and sources. The same layer is available to the final sequence report.

Current curated entry: **Gallus gallus / GeneID 101748883 (LOC101748883)** is
displayed as **NANOG-like**, retaining the original LOC identifier and the
historical names Nanogl1 / Nanog-like 1. This is an exact species/source-gene rule;
it is not applied to other species by sequence similarity. See the
[source gene](https://www.ncbi.nlm.nih.gov/gene/101748883) and the original
[description of the adjacent chicken Nanogl1 locus](https://doi.org/10.1002/dvdy.20927).
The label does not establish expression, protein production or a pluripotency
function, and does not classify the locus as an inactive pseudogene.

Other named annotated genes retain their provided names even if database hits
have different names. Unnamed LOCs and predictions supported only by homology
are explicitly **family candidates**, with copy identity unresolved. These
names are not propagated to other species or used to infer copy orthology.

A prediction inherits an annotation label only if genome, contig, strand, full
phased CDS blocks, stop blocks and an available CDS sequence fingerprint agree
exactly with one unambiguous annotated gene. Overlap, equal proteins, physical
group membership and a common best hit are insufficient.

Saved copy labels and comparison groups remain user decisions. Distinct copies
can intentionally share a family-level comparison group for all-by-all sequence
comparison; the group name must not be presented as each copy's specific identity.
No selection, sequence, alignment or evolutionary assignment is changed by a
presentation refresh. Future curated entries require explicit source-gene and
species provenance, rather than gene-name substring matching.
