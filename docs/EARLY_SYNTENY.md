# Automatic early context evidence

The early-synteny backend prepares annotation-supported context evidence from the query, study-genome list, and supplied source datasets. It does not ask the reviewer to assemble neighbour panels, translate CDS, create indexes, or match protein IDs manually. The normal manual review stop remains a separate decision after computation.

## Inputs

A dataset TSV uses `genome_id` and optional `gff3`, `protein`, `cds`, `species`, `assembly_accession`, `annotation_provenance`. Paths are relative to that TSV. Supply GFF3 plus protein/CDS FASTA, or GFF3 plus genomic `fasta` for protein derivation. With supplied protein/CDS and GFF sequence-region lengths, the genomic FASTA is neither required nor read/hashed; it can still be downloading. Missing contig lengths remain explicit. An explicit source-dataset TSV can add reference assemblies absent from the study panel; otherwise supplied study datasets also serve as possible sources. A source/target row with the same genome ID may fill missing fields of the other row, but conflicting nonempty paths are rejected.

When explicit annotation/protein/CDS paths are absent, only exact NCBI sibling names (`*_genomic.gff[.gz]`, `*_protein.faa[.gz]`, `*_cds_from_genomic.fna[.gz]`) are discovered. Other annotation releases are not searched recursively or downloaded. Missing annotation is reported as `annotation_unavailable` and is never interpreted as absence of the target gene.

Source resolution checks versioned protein/transcript identifiers (including those in FASTA IDs), supplied GeneIDs, and exact product-sequence hashes. All usable annotation products participate in this binding. A supplied GeneID must also have a matching product sequence. With no identifier match, an exact sequence can identify possible source loci. Multiple matches remain explicitly multiple; an identifier/sequence conflict is not silently resolved by a similar sequence elsewhere. Isoforms of one resolved source gene share one reference context. Query members not represented by the provided source datasets remain `source_locus_unresolved` or `identifier_sequence_mismatch`; other references can proceed.

Protein queries and translated spliced CDS queries are supported. An unspliced DNA query does not receive fabricated protein-reference provenance and is reported with `coding_reference_unavailable`.

## Automatic preparation

1. Parse each supplied GFF3 once into distinct gene positions and parented CDS models.
2. Map an available protein FASTA or spliced CDS FASTA to exact GFF3 accessions. NCBI CDS headers are additionally checked against contig, GeneID and native CDS intervals when available: the same protein accession can occur on primary/alternate loci. Unresolved contradictory CDS sequences quarantine only the affected transcript in `unavailable_products.json`; they do not abort an otherwise usable proteome or select an arbitrary sequence. CDS translation uses the standard nuclear code; explicitly different translation tables are excluded. If neither is supplied, derive proteins from GFF3 CDS and the source genomic FASTA using the existing coding-sequence checks.
3. Keep all usable products for source resolution, and one longest usable product per gene for the full target search catalogue. Missing/unusable proteins remain visible as incomplete evidence.
4. Export every usable product as `all_products.faa`, searching identical sequences once, with a binding row for every transcript/gene in `products.jsonl` and `products.tsv`. Protein-first discovery must use these files, not the longest-product catalogue. Missing/unusable models are listed in `unavailable_products.json`.
5. Retrieve five coding neighbours per genomic side, excluding genes overlapping the target interval and other resolved source target genes. Do not count isoforms as independent neighbours.
6. Reserve neighbours up to twenty per side within a one-megabase extent, and search their unique proteins in one BLASTP panel against each complete available target catalogue. No target-gene BLASTP, target MSA, genome-wide neighbour TBLASTn, or Tiberius is run here.
7. Reuse the resulting protein correspondences for every candidate locus. Expand annotation lookups from five to ten/twenty when a two-sided ordered context has not been found; this does not rerun BLAST.

BLASTP runs against the full available protein catalogue, not just candidate-window proteins. Default preliminary-match gates are 50% coverage of each protein and 30% identity. All passing hits are examined; hits scoring within 95% of the global best remain ambiguous. Several source genes matching one target gene are excluded from independent-anchor counts. These settings are inspection aids, not an orthology classifier or universal biological thresholds.

Order is evaluated in each target gene's transcriptional orientation. Matched upstream/downstream sides and relative gene orientations are retained separately. Two-sided support requires at least two unique preliminary anchors, at least one on each matching side, and consistent order/orientation. Rearrangements, missing annotations and incomplete assemblies can prevent support for a true locus. No candidate is removed, approved, or relabelled automatically.

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

Each reference retains `reference_id`, source genome/gene IDs, original query aliases, source-resolution methods, coordinates, strand, neighbours, limits and actual window. Neighbours have gene/protein IDs, coordinates, genomic side, target-relative side/orientation, nearest-gene rank and distance.

Each locus retains the original `genome_id` and `locus_id`, `status`, requested/used neighbour count, actual window, `neighbors`, `comparisons`, limits and flags. Comparisons contain `reference_id`, `shared_anchors`, `left_anchors`, `right_anchors`, `order_concordance`, `orientation_concordance`, `matches`, `ambiguous_matches`, `colliding_matches`, and flags. A match contains source/target gene IDs and names, sides, starts, relative strands, protein IDs, bit score, E-value, identity and both coverages.

Locus statuses include `two_sided_context_support`, `one_sided_context_support`, `order_or_orientation_differs`, `insufficient_anchors`, `source_locus_unresolved`, `annotation_unavailable`, `annotation_contig_missing`, `annotation_contig_length_mismatch`, `protein_matching_unavailable`, and `dataset_error`. Match support is always preliminary protein-homology context evidence, not established orthology.

## All-product catalogue API for protein-first discovery

`build_catalog(row, cache_base, runner, known_fasta_digest=None)` returns `all_products_fasta`, `products_tsv`, `products_jsonl`, `all_product_models`, `unique_product_sequences`, `unavailable_product_models`, `cache`, `key`, `genes`, `contigs`, `inputs`, and preparation/reuse metadata. `iter_products(catalog)` streams binding records. Existing `proteins.faa` retains one longest usable product per gene solely for neighbour matching.

Each product binding has `search_id`, `sequence_sha256`, `protein_ids` (original accessions), `model_id`, `transcript_ids`, `gene_id`, `ncbi_gene_id`, `gene_name`, `contig`, `gene_start`, `gene_end`, CDS `start`/`end`, `strand`, `protein_aa`, native `cds` blocks with phase, `stop_codons`, `origin`, annotation `exception`, and `aliases`. Coordinates are zero-based half-open. Several products/genes can share a search ID when the complete product sequences are identical; no binding is removed. Protein BLAST coordinates remain amino-acid coordinates and must not be presented as genomic alignments.

## Compact review display

The default diagram shows the candidate once and every available reference locus, with five extracted neighbors per side and equal schematic spacing. Variants from one annotated source locus share one reference row. All rows are oriented along their own target gene; matching colors propose a link to the same candidate neighbor without asserting orthogroups. Same-assembly/self-control rows are labeled and do not receive independent-support colors. Native coordinates and linked query IDs remain in tooltips. The optional distance scale, ±50 kb view and single-reference inspection live under **How this is checked, missing references and pair details**. Changing display scope only exposes saved neighbors; it does not search or predict additional genes. A per-query table explains missing reference neighborhoods. The backend already calculates comparisons against all resolved reference loci, independent of what the UI displays.

Dashed links in optional single-reference inspection are preliminary protein matches. A different symbol (for example FOXJ2 versus FOXJ3) carries a gene-family/paralog review note, while matching symbols still do not establish orthology. Self-controls and other same-assembly comparisons are labeled and excluded from cross-species summary counts. Legacy heuristic status codes remain in technical details for reproducibility.

The reviewer may shortlist proposed pairs or individual reference/candidate genes, including genes without a match. Draft JSON preserves these choices and `neighbors.to_inspect.tsv` exports exact locus/reference/gene IDs. The action is `inspect`, not approval. Target-window export remains separate; downstream anchor validation does not yet consume this shortlist automatically. See [metric definitions](STAGE1_METRICS.md).
