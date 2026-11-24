# Automatic model placement and Stage4 review

The same gene-independent rule assesses every exact model on its shared gene
tree. There is no scoring leaf picker, branch picker, or editable T value.
Sequence membership and the biological outgroup remain scientific choices at
the tree-building step. The reviewer can inspect and change the outgroup at
Stage4. These choices never select a convenient branch support for a score.

## Calculation

An independent species Newick or a root-to-tip taxonomic lineage table supplies
the expected relationships. The bundled 18-species classification scaffold is
an explicit study input, independent of the gene tree. It retains polytomies,
including unresolved within-order relationships; it is not inferred from the
protein sequences and is not presented as a universally resolved species tree.
Any other taxon set can use --species-tree or --species-lineages. No gene names
or vertebrate-specific cases occur in the calculation.

For each focal model, enumerate triples of other species with one unambiguous
eligible primary protein per species. Compare the resulting unrooted quartet
with the independent scaffold. A species with multiple distinct eligible peer
proteins is ambiguous and is not silently assigned its best-fitting copy.
Alternative focal models can each be assessed, replacing their own species'
primary in the comparison. References and isoforms are not independent species.

For an assessed quartet, use the minimum reported support on its induced central
path. A concordant quartet contributes support/100; a conflicting quartet
contributes zero. The proposed component is:

    T = 15 × sum(concordant quartet contributions) / assessed quartet count

Only resolved, distinct-tip quartets with known support are assessed. Unknown
species, absent focal models, identical-sequence ties, unresolved relationships
and unknown support fields are reported separately. They are not measured
zeros. No assessed quartet means NA. Coverage, concordance, discordance, mean
support and informative counts accompany every proposal. A partial estimate is
labelled partial. For panels exceeding 2000 possible quartets per focal species,
a deterministic taxon-based sample of 2000 is used and declared explicitly;
ordinary 18-species panels use exhaustive enumeration.

Newick root location, sibling order, gene names, model identifiers and display
outgroup do not determine the score. Multiple support fields are never guessed:
provide --support-field when the input uses a declared multi-support format.
This metric is a versioned descriptive rule, not IQ-TREE's gCF or sCF, not a
calibrated orthology probability, and not a diagnosis of annotation error.
Duplication, loss, incomplete lineage sorting and genuine gene history can
produce discordance. Quartets share sequences and are not independent tests.

## Review gate

A calculated proposal is visible before approval but is excluded from the total.
The Stage4 reviewer checks the exact tree, sequence set and outgroup (or knowingly
accepts an unrooted interpretation), then confirms the tree. The confirmation
admits the automatically calculated components for that tree; there is no manual
selection of a convenient split. A reviewed model can still have T=NA.

Confirmation is bound to tree and alignment identity, membership, outgroup and
its basis, species scaffold, protein identities, method and code version. Any
change invalidates it. Reopening or exporting a working session preserves the
bound decision. An apply operation validates the returned data and recomputes
proposals; client-supplied scores do not override the calculation.

## Commands

    gene-context tree-assess --tree gene.nwk --models models.json \
      --species-tree species.nwk --out proposed

models.json is a list of records containing model_id, tree_id and taxon; optional
protein_sha256 and primary fields preserve identities and peer eligibility.
The output is a proposal, not an automatically approved study score.

    gene-context stage4 --handoff Stage4-input.json --out new-dossier \
      --species-tree species.nwk --structure-provider none
    gene-context tree-review --results new-dossier
    gene-context tree-review-apply --results new-dossier \
      --reviews returned-tree-reviews.json --out confirmed-dossier

An older sealed native Stage4 dossier can be prepared into a separate review copy
with tree-review --results OLD --out NEW. Original scientific inputs are retained.
A prepared review bundle can also launch with python3 START_REVIEW.py after
extraction; no research server or local Go compiler is required.

Method background: [IQ-TREE branch support](https://www.iqtree.org/doc/Frequently-Asked-Questions#how-do-i-interpret-ultrafast-bootstrap-ufboot-support-values),
[concordance factors](https://www.iqtree.org/doc/Concordance-Factor), and
[NCBI taxonomy as classification](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/).
The proposed formula above is this tool's own documented rule.
