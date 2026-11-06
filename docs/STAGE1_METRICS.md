# Reading stage 1 evidence

The compact view shows one ranked candidate and its neighbor diagram. More candidates appear ten at a time. The sequence metrics describe the **best-scoring single query/product chain**; open **Sequence matches and metric explanations** to compare all supplied queries and annotated isoforms. Click a metric card or focus a table heading marked ⓘ for its definition.

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

- Purple: target gene. Both target genes point right; neighbor arrows show direction relative to it.
- Gray: displayed neighbor without a retained proposed pair in both visible rows.
- Colored arrows and dashed links: **proposed protein correspondence, orthology unverified**. Matching names alone do not validate a pair; different names flag possible gene-family/paralog issues.
- Default: five extracted neighbors per side, equal spacing. This is gene order, not a physical DNA scale. Native coordinates and distances are retained in tooltips and detailed tables.
- **Change reference or view**: select the reference, genomic-distance scale, ±50 kb window or larger saved neighborhood. The ±50 kb view is measured from each target's edges and contains only the already-extracted neighbor inventory.
- The caption reports pairs visible in both rows versus all proposed pairs. A pair outside the view is not missing from the evidence.
- Comparisons with the same locus are self-controls. Other same-assembly comparisons are labeled too. Neither contributes to the cross-species overview.

## Detailed neighbor metrics

These measures are collapsed under **Inspect proposed pairs → Matching rules, ambiguity and legacy diagnostic measures**. They describe the proposed correspondences, which may contain mistaken matches.

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
