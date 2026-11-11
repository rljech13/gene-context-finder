# Automatic early context evidence

The early-synteny backend prepares annotation-supported context evidence from the query, study-genome list, and supplied source datasets. It does not ask the reviewer to assemble neighbour panels, translate CDS, create indexes, or match protein IDs manually. The normal manual review stop remains a separate decision after computation.

## Inputs

A dataset TSV uses `genome_id` and optional `gff3`, `protein`, `cds`, `species`, `assembly_accession`, `annotation_provenance`. Paths are relative to that TSV. Supply GFF3 plus protein/CDS FASTA, or GFF3 plus genomic `fasta` for protein derivation. With supplied protein/CDS and GFF sequence-region lengths, the genomic FASTA is neither required nor read/hashed; it can still be downloading. Missing contig lengths remain explicit. An explicit source-dataset TSV can add reference assemblies absent from the study panel; otherwise supplied study datasets also serve as possible sources. A source/target row with the same genome ID may fill missing fields of the other row, but conflicting nonempty paths are rejected.

When explicit annotation/protein/CDS paths are absent, only exact NCBI sibling names (`*_genomic.gff[.gz]`, `*_protein.faa[.gz]`, `*_cds_from_genomic.fna[.gz]`) are discovered. Other annotation releases are not searched recursively or downloaded. Missing annotation is reported as `annotation_unavailable` and is never interpreted as absence of the target gene.

Source resolution checks versioned protein/transcript identifiers (including FASTA IDs), supplied numeric GeneIDs, and exact product-sequence hashes. All usable annotation products participate. Exact identifier/sequence bindings are preferred. A supplied curated **GeneID** may also bind a source annotation locus for neighborhood lookup even when the query isoform is absent or differs in this annotation release; this has status `gene_id_linked_sequence_differs`, `context_binding_only=true`, and `query_sequence_verified=false`. The discrepancy remains recorded in `mismatching_sources`. This binds the gene's context; it does not verify the query sequence, replace it with an annotated product, or approve any candidate. Gene symbols alone never bind source loci.

Multiple annotation locations remain explicitly multiple. Without a curated GeneID, identifier/sequence conflicts remain unresolved; a similar protein is not used to invent a source. With no identifier conflict, an exact sequence can identify source loci. Isoforms of one resolved source gene share one reference context. Missing source datasets remain `source_locus_unresolved`; source datasets can be supplied independently of target genomes using `--source-datasets`.

Protein queries and translated spliced CDS queries are supported. An unspliced DNA query does not receive fabricated protein-reference provenance and is reported with `coding_reference_unavailable`.

## Automatic preparation

1. Parse each supplied GFF3 once into distinct gene positions and parented CDS models. Alignment-evidence rows (`cDNA_match`, `match_part` and related match types) are ignored; their legal multipart IDs are not interpreted as duplicated gene/transcript identifiers. Duplicate annotation gene/transcript IDs remain errors.
2. Map an available protein FASTA or spliced CDS FASTA to exact GFF3 accessions. NCBI CDS headers are additionally checked against contig, GeneID and native CDS intervals when available: the same protein accession can occur on primary/alternate loci. Unresolved contradictory CDS sequences quarantine only the affected transcript in `unavailable_products.json`; they do not abort an otherwise usable proteome or select an arbitrary sequence. CDS translation uses the standard nuclear code; explicitly different translation tables are excluded. If neither is supplied, derive proteins from GFF3 CDS and the source genomic FASTA using the existing coding-sequence checks.
3. Keep all usable products for source resolution, and one longest usable product per gene for the full target search catalogue. Missing/unusable proteins remain visible as incomplete evidence.
4. Export every usable product as `all_products.faa`, searching identical sequences once, with a binding row for every transcript/gene in `products.jsonl` and `products.tsv`. Protein-first discovery must use these files, not the longest-product catalogue. Missing/unusable models are listed in `unavailable_products.json`.
5. Retrieve five coding neighbours per genomic side, excluding genes overlapping the target interval and other resolved source target genes. Do not count isoforms as independent neighbours.
6. Reserve neighbours up to ten per side within a one-megabase extent, and search their unique proteins in one BLASTP panel against each complete available target catalogue. No target-gene BLASTP, target MSA, genome-wide neighbour TBLASTn, or Tiberius is run here.
7. Reuse the resulting protein correspondences for every candidate locus. Expand annotation lookups from five to ten per side when a cross-species two-sided ordered context has not been found; this does not rerun BLAST.

BLASTP runs against the full available protein catalogue, not just candidate-window proteins. Default preliminary-match gates are 50% coverage of each protein and 30% identity. All passing hits are examined; hits scoring within 95% of the global best remain ambiguous. Several source genes matching one target gene are excluded from independent-anchor counts. These settings are inspection aids, not an orthology classifier or universal biological thresholds.

Order is evaluated in each target gene's transcriptional orientation. Matched upstream/downstream sides and relative gene orientations are retained separately. Two-sided support requires at least two unique preliminary anchors, at least one on each matching side, and consistent order/orientation. Rearrangements, missing annotations and incomplete assemblies can prevent support for a true locus. The maximum is ten genes per side; larger values are rejected. Same-assembly comparisons are retained as controls and excluded from the expansion stopping rule and locus-level cross-species support status. With only controls, the status is `within_assembly_controls_only`. No candidate is removed, approved, or relabelled automatically.

## Cache and provenance

The caller supplies a writable cache outside the checkpoint (on `/mnt/data` for host2). GFF3, FASTA and supplied protein/CDS content hashes, backend source hash, tool versions and matching parameters identify reusable artefacts. OS-released file locks protect concurrent writers. New entries are published atomically; reuse checks their recorded output hashes. Input byte counts/mtime are checked for changes during construction. Annotation protein databases are cached independently of a particular query panel.

An optional `annotation_provenance` JSON sidecar can record a compact annotation snapshot. Its genome ID and snapshot SHA256 must match the actual GFF, and its original source SHA256 must match the recorded expected fingerprint. That derivation is embedded in the dataset receipt. Compact pilot snapshots preserve sequence-region directives and all gene/mRNA/CDS/stop-codon coordinates, phases and parent/accession fields; descriptive and exon/UTR rows may be omitted.

Exact input hashing still requires input reads. Cache reuse avoids GFF parsing, protein extraction, database construction and repeated alignment, but is not claimed to eliminate all disk I/O.

The output is a separate evidence directory containing `evidence.json`, `summary.tsv`, `reference_neighbors.faa`, and per-target raw matching evidence under `matches/`. A caller must include that directory in its checkpoint seal. A portable review does not need the cache to display these saved results.

## Integration API

```python
from genecontext.early_synteny import run_early_synteny

payload = run_early_synteny(
    data,                              # existing stage1 evidence dictionary
    genomes_manifest=study_datasets,
    source_genomes_manifest=source_datasets,  # optional
    cache_dir=shared_cache,
    out_dir=output / 'raw' / 'early_synteny',  # must not exist yet
    threads=8,
)
```

Run this once in the parent after discovery workers merge, before the final report and seal. It does not mutate `data` or its source checkpoint. Attach `payload['loci']` by the exact `(genome_id, locus_id)` pair. The global `references`, `source_resolution`, `datasets`, `errors`, parameters and command receipts remain available to the report. A future input file selection may change paths, but does not require manually constructing a reference panel.

A standalone pilot can enrich saved discovery without another genome search:

```bash
python -m genecontext.early_synteny --results SAVED_STAGE1 \
  --genomes DATASETS_TSV --source-genomes SOURCE_DATASETS_TSV \
  --cache SHARED_CACHE --out NEW_OUTPUT --threads 8
```

On host2 use the existing constrained launcher and explicit source snapshot/PYTHONPATH. This command never changes `current`, existing review IDs, evidence, or decisions.

## Output schema

Top level: `schema=early_synteny_v1`, `status`, `parameters`, `source_resolution`, `datasets`, `references`, `loci`, `errors`, `commands`, `elapsed_seconds`, `tools`.

Each reference retains `reference_id`, source genome/gene IDs, original query aliases, source-resolution methods, coordinates, strand, neighbours, optional display-only annotations, limits and actual window. Neighbours have gene/protein IDs, coordinates, genomic side, target-relative side/orientation, nearest-gene rank and distance.

Each locus retains the original `genome_id` and `locus_id`, `status`, requested/used neighbour count, actual window, `neighbors`, `display_only_genes`, `comparisons`, limits and flags. Comparisons contain `reference_id`, `shared_anchors`, `left_anchors`, `right_anchors`, `order_concordance`, `orientation_concordance`, `matches`, `ambiguous_matches`, `colliding_matches`, and flags. A match contains source/target gene IDs and names, sides, starts, relative strands, protein IDs, bit score, E-value, identity and both coverages.

Locus statuses include `two_sided_context_support`, `one_sided_context_support`, `order_or_orientation_differs`, `insufficient_anchors`, `source_locus_unresolved`, `annotation_unavailable`, `annotation_contig_missing`, `annotation_contig_length_mismatch`, `protein_matching_unavailable`, and `dataset_error`. Match support is always preliminary protein-homology context evidence, not established orthology.

## All-product catalogue API for protein-first discovery

`build_catalog(row, cache_base, runner, known_fasta_digest=None)` returns `all_products_fasta`, `products_tsv`, `products_jsonl`, `all_product_models`, `unique_product_sequences`, `unavailable_product_models`, `cache`, `key`, `genes`, `contigs`, `inputs`, and preparation/reuse metadata. `iter_products(catalog)` streams binding records. Existing `proteins.faa` retains one longest usable product per gene solely for neighbour matching.

Each product binding has `search_id`, `sequence_sha256`, `protein_ids` (original accessions), `model_id`, `transcript_ids`, `gene_id`, `ncbi_gene_id`, `gene_name`, `contig`, `gene_start`, `gene_end`, CDS `start`/`end`, `strand`, `protein_aa`, native `cds` blocks with phase, `stop_codons`, `origin`, annotation `exception`, and `aliases`. Coordinates are zero-based half-open. Several products/genes can share a search ID when the complete product sequences are identical; no binding is removed. Protein BLAST coordinates remain amino-acid coordinates and must not be presented as genomic alignments.

## Compact review display

The report uses one all-reference diagram, with up to ten saved independent neighbors per side. Scope/scale controls are beside it. Equal spacing shows order; the distance scale shows kilobases relative to the oriented target. The actual displayed span is separate from the later ±50 kb window. Historical results containing more neighbors are retained in evidence but capped in the view.

`display_only_genes` retains other resolved query-source genes and overlapping coding annotations excluded from neighbor anchors. These are visible with purple outlines, consume no independent-neighbor slots, and never enter the neighbor protein panel or support/expansion calculations. The central anchor is not duplicated. Historical reports can recover known omitted source copies from their saved reference list, without recomputation. Display-only genes are not a claim that an old annotation neighborhood was exhaustively reconstructed.

Solid curves link saved protein proposals through the same candidate neighbor. Dashed curves include same-assembly controls; amber dotted curves include unresolved assignments. No curve bridges an intervening row without that correspondence, and references have not been searched directly against one another. Exact same-assembly annotation identities can be shown as controls, but display-only target copies are excluded. Gray annotations remain visible; lack of a link does not establish absence. Hovering traces a candidate-linked chain.

The optional match table retains protein identity, both coverages, bit score, E-value and ambiguity. Neighbor/pair follow-up requests persist in draft JSON and export separately; `stage1-select` still consumes target windows. Source availability is listed once in the global query-source table. There is no second pair diagram or genome-wide best-neighbor metric attached to an unrelated selected locus. See [metric definitions](STAGE1_METRICS.md).
