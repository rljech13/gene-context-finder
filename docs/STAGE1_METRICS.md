# Reading stage 1 evidence

The review unit is a **candidate genomic locus**. Its compact card shows a row for every supplied query variant: query/source, covered positions and percentage, aligned identity, and chain score. A row summarizes that query's best-scoring single query/product chain at this locus. Weak matches stay visible and are marked **Below filters** when no non-ambiguous chain from that query passes the current thresholds. A query without a retained chain is explicitly marked **No retained match**. Identical input sequences retain separate variant names but share evidence; the variant count is not a count of independent votes.

More candidate loci appear ten at a time. Open **Sequence matches and metric explanations** for all annotated isoforms and HSPs; the expanded metric cards explicitly describe the single strongest query. Click a metric card or focus a table heading marked ⓘ for its definition.

| Metric | Meaning | What it does not establish |
|---|---|---|
| Query coverage | Union of aligned query positions / full query length. Overlapping positions count once. | Completeness of the candidate gene, or identity of unaligned sequence. |
| Aligned identity | Estimated fraction of identical residues in retained local HSPs. Target-search query overlaps are prorated. | Identity across the entire query. |
| Subject coverage | Union of aligned positions / annotated candidate protein length, for BLASTP. | Completeness of the supplied annotation. |
| Chain score | Sum of HSP bit scores prorated by fresh query positions. Ranking uses the strongest single-query chain. Scores from query variants or isoforms are not added. | Probability, orthology or synteny. |
| Best HSP E-value | Smallest local HSP E-value, not a combined gene-level statistic. Smaller means less expected chance similarity in that search. A reported 0 is numerical underflow. | Probability that this is the correct ortholog. |
| Rank | Order within one genome/search group by chain score. | Automatic gene assignment. |
| Display gate | One non-ambiguous chain passes coverage, identity and E-value gates together. It may differ from the summary's strongest chain. | Biological approval or evidence that a filtered candidate is absent. |

For example, the saved chicken NANOG candidate covers **299 / 449 aa = 66.6%** of the strongest query, at **100% identity in the aligned parts**. These values are compatible: the remaining query positions are not part of the identity denominator.

The meaning of local alignment identity, bit scores and E-values follows the [NCBI BLAST glossary](https://www.ncbi.nlm.nih.gov/books/NBK62051/) and [NCBI explanation of E-values](https://www.nlm.nih.gov/ncbi/workshops/2023-08_BLAST_evol/e_value.html). The chain aggregation and ranking rules above are specific to this pipeline.

## Neighbor diagram

The default overview contains the candidate locus once and **all available query-source loci**, one row each. Query variants from the same source gene share a row. Every row preserves its own gene order, with its target pointing right.

- Purple: target gene. Arrows show direction relative to it.
- Matching colors: proposed protein correspondence to a particular candidate neighbor, **orthology unverified**. Colors are local diagram keys, not confirmed orthogroups. A candidate gene is colored only if a corresponding reference gene is also visible. Gray does not mean the gene is biologically absent.
- Default: five saved neighbors per side, equal spacing. Native coordinates, distances and linked query IDs are available in tooltips. The count in each reference row distinguishes visible from all retained proposed pairs.
- Self-control and other same-assembly rows are labeled and do not receive independent-support colors. They are excluded from the cross-species overview count.
- **How this is checked, missing references and pair details** contains the method, per-query source availability, optional genomic-distance/±50 kb/larger-neighborhood controls and single-reference inspection. Changing the view uses saved evidence only.
- The detailed two-row inspection retains dashed links, individual protein metrics, and neighbor-selection controls. This display can focus on one reference without hiding other references from the main overview.

A FASTA query has no genomic neighborhood by itself. Its source locus must be resolved in supplied annotation/coding-product datasets. Missing source datasets, identifier/sequence conflicts and unresolved source loci are shown explicitly for each query. An available source neighborhood can be shared by multiple isoforms; neither isoform counts nor reference counts are independent votes for orthology.

## Detailed neighbor metrics

These measures are collapsed under **How this is checked, missing references and pair details → Inspect one reference → Inspect proposed pairs → Matching rules, ambiguity and legacy diagnostic measures**. They describe the proposed correspondences, which may contain mistaken matches.

| Metric | Numerator / denominator or rule |
|---|---|
| Proposed pairs | Retained non-competing, non-colliding protein correspondences in the extracted neighborhoods. Default gates: ≥50% coverage of each protein, ≥30% identity, E ≤ 1e-5. These gates are independent of target display filters. |
| Neighbor coverage / identity / score | Calculated from non-overlapping BLASTP HSPs; an HSP overlapping either an existing query or subject span is omitted. Identity is weighted by aligned length; score sums the retained bit scores. One longest usable protein per neighboring gene is used. |
| Same target-relative side | Pairs upstream in both rows plus pairs downstream in both rows / all proposed pairs. Undefined with no pairs. |
| Relative strand agreement | Pairs with the same direction relative to their target / all proposed pairs. Undefined with no pairs. |
| Pairwise order agreement | Concordant gene-start orders / all pairs of proposed correspondences, after orienting both target genes to the right. For n correspondences the denominator is n(n−1)/2; undefined for n < 2. |
| Competing assignments | More than one target reaches 95% of the best score in the full available target protein catalogue. Retained separately, excluded from proposed-pair counts. |
| Colliding assignments | Several reference neighbors propose the same candidate gene. Excluded from independent pair counts. |

No single count or legacy status code is presented as a curated synteny conclusion. Missing annotations, unusable products and limited extracted neighborhoods can all leave gaps.

## Neighbor choices

The optional checkboxes request further inspection of proposed pairs or individual genes, including unmatched neighbors. Draft JSON preserves these choices; **Download neighbor shortlist TSV** exports `neighbors.to_inspect.tsv` with review, genome, locus, reference and gene IDs. This does not approve orthology. `stage1-select` currently exports selected **target-locus windows** only; downstream anchor validation is not yet connected to this shortlist.
