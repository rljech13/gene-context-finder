# Stage 1: one gene, multiple reference variants, candidate loci

## Contract

Input: a manually selected reference FASTA for **one gene/search group**, its metadata and an explicit list of genomes. Include all chosen annotation variants; the reference set need not contain exactly four sequences. No genome annotation is required.

Output: candidate genomic loci with all supporting query alignments, a manual decision file, and explicitly selected intervals for subsequent synteny analysis. This stage does not determine orthology, annotate neighbouring genes or measure final cross-species conservation.

The new `stage1` commands implement this contract. The older `search` command has a different, model-first checkpoint format. They must not be mixed. `stage1-select` currently exports TSV/BED regions; these exports are **not yet accepted directly by the legacy `contexts` command**. The next-stage interface should consume the approved intervals without searching for targets again.

## First test

Run from the repository root with BLAST+, samtools and miniprot on PATH:

```bash
gene-context stage1 --gene POU5 \
  --query examples/queries/pou5_test/proteins.faa --query-type protein \
  --query-metadata examples/queries/pou5_test/queries.tsv \
  --genomes examples/queries/pou5_test/controls.genomes.tsv \
  --refine miniprot --out results/pou5_stage1 --threads 4
```

Open `results/pou5_stage1/report.html`. This test searches **small annotated reference windows**, not 18 complete genomes. For a complete run, replace only `--genomes` with the TSV of the actual assemblies. Include source-reference assemblies as ordinary genome rows if source controls should run alongside the study genomes.

On host2, submit this command through the constrained launcher described in [LIMITED_EXPERIMENTS.md](LIMITED_EXPERIMENTS.md), using `/usr/bin/python3 -m genecontext` and absolute `/mnt/data` paths. The deployment must contain this version; a v0.5.2 source snapshot does not provide `stage1`. Do not bypass the shared 250 GiB RAM cap.

## Inputs

Genome TSV, with paths relative to that TSV or absolute:

```tsv
genome_id	fasta
human	assemblies/human_genomic.fna.gz
chicken	assemblies/chicken_genomic.fna.gz
```

One row means one assembly. The stage searches every listed assembly, including source species; there is no self-hit exclusion. It does not recursively treat CDS/transcript FASTAs as additional genomes. Gzip is accepted. Genome errors remain in the report and cause a nonzero final exit code.

Query metadata must assign every FASTA ID exactly once. Minimum useful columns are `query` and `species`; optional columns include `group`, `gene_symbol`, `transcript` and source links. If `group` is supplied, it must match `--gene`. Without species metadata, species support is unknown (reported as zero labelled species), not inferred from filenames.

Identical sequences are searched once while preserving every transcript/species alias. Variant count, distinct sequence count and supporting species count are separate quantities. A larger number of annotated isoforms does not multiply a locus's ranking score.

`--query-type protein` uses TBLASTn. `cds` translates the supplied spliced CDS and uses TBLASTn; this is **not an additional BLASTN run**. `dna` uses BLASTN and does not perform protein-guided reconstruction. Protein/CDS validation assumes the standard nuclear genetic code and the documented query alphabet. Use separate runs for protein and nucleotide evidence.

## Automatic assembly before manual review

1. BLAST retrieves local matches with the relatively permissive `--search-evalue` ceiling. Raw output is retained. All target contigs can be reported; no default top-500-contig cap is used.
2. Exact duplicate HSPs are removed. The seed length and identity gates are applied; rejected seeds and reasons remain in raw JSON.
3. Compatible HSPs are chained on the same contig and strand in query/transcriptional order. Query-coordinate resets prevent a simple join of tandem full-length copies. Intron length, overlap tolerance and total locus span are configurable. The algorithm extracts best collinear paths and retains residual chains; it does not enumerate every possible exon chain.
4. Chains from different queries are initially grouped by shared **aligned genomic blocks**, not merely overlapping envelopes or distance.
5. Local miniprot reconstruction tests each proposal. The full GFF3 is retained; the best coverage-times-identity model per unique query is offered in the review. The mapper's alternative-output limit is explicit, default 10 per query/local window; it does not limit the number of BLAST loci.
6. Spliced models reconcile fragmented proposals into loci. A fixed representative provides the proposed coding anchor; weaker chains never extend it across a neighbouring copy. A chain supported by separate model anchors is split by its HSP-to-anchor assignments. All retained HSPs stay in the evidence. Ambiguous assignments are flagged and do not satisfy review thresholds by themselves. Unassigned fragments remain BLAST-only candidates.

The initial and model-based grouping rules differ intentionally: model blocks help resolve fragments that BLAST extends into intronic sequence. Model consolidation uses fixed representatives, avoiding transitive extension through bridging matches. Nearby distinct genes can remain distinct.

This is a heuristic for candidate loci, not a validated gene caller. A wrong or missing spliced model can still leave a split/merged proposal. Miniprot failures preserve BLAST candidates. `--refine miniprot` requires the executable; default `auto` records its absence and leaves BLAST-only proposals. `--refine none` skips this step explicitly.

The [host2 validation receipt](benchmarks/stage1_pou5_v060_2026-11-06.json) records the reference-window outcomes and tests. It is not a complete-genome sensitivity benchmark.

The default anchor is the representative's predicted **coding interval**, not an established full gene/UTR interval. Without reconstruction it is the BLAST-chain envelope. Both can be corrected manually.

## Review and filters

The HTML report works offline. It has text search, coverage/identity/E-value display filters, compact alignment tracks and a switch to show every query variant. The reference-support table preserves transcript identities. Alternative source-derived CDS models can be used to set anchor bounds.

Default review thresholds are `--min-query-coverage 0.15`, `--min-identity 0.25`, `--max-evalue 1e-5`. **One and the same reference chain must pass all three.** Identity is measured within aligned sequence; query coverage uses the union of covered query positions, never the sum of overlapping matches. For small overlaps, identity and score contributions are prorated by the fraction of newly covered query positions; aggregate identity is an estimate, not a new residue-level alignment. These are BLAST-chain metrics, distinct from the local model's alignment metrics.

No candidate is deleted by these review thresholds. `loci.tsv` contains all retained loci; `passing_loci.tsv` is the saved passing subset. Neither constitutes approval. Browser filtering changes the view only.

To save different review thresholds without rerunning BLAST, chaining or reconstruction:

```bash
gene-context stage1-filter --results results/pou5_stage1 \
  --min-query-coverage 0.50 --min-identity 0.35 \
  --out results/pou5_stage1_stricter
```

This creates a new review ID and blank decisions. It does not transfer old approvals. Review E-value cannot exceed the original search ceiling.

Search/assembly settings require a fresh `stage1` run: `--search-evalue 1e-3`, `--seed-min-length 20` (aa for protein/CDS, nt for DNA), `--seed-min-identity 0.20`, `--max-intron 200000`, `--max-locus-span 2000000`, `--chain-max-overlap 0.15`, `--locus-shared-fraction 0.5`, `--refine-padding 50000`, `--refine-max-models 10`. These are pilot defaults, not universal biological thresholds. A raw match omitted by BLAST or the seed gates cannot be recovered by display filtering.

## Explicit decisions and synteny intervals

In the report, record `inspect`, `exclude` or `unresolved` for **every locus**. Genomes without candidates also have a required decision row. Notes are required for exclusions/unresolved cases, manual boundary/window changes and below-filter selections. A failed genome cannot be marked `inspect` without a locus.

Download `loci.decisions.tsv` from the report. Its export contains **all rows, including hidden loci**. A saved decision file can be reloaded in the same report. Browser-local storage is a convenience; the downloaded TSV is the portable review record.

`anchor_intervals` uses zero-based `start:end`. Change it to correct boundaries. To split an erroneously joined proposal, use `start1:end1;start2:end2` and explain the change in notes. Each anchor must overlap its original candidate and stay within the contig. Several independently discovered loci can be selected. A merged inspection window can be specified explicitly using `window_start` and `window_end` on the chosen locus, with other rows still receiving explicit decisions.

By default, the proposed synteny window is anchor ±50 kb, clipped to the contig. A custom window must contain the anchor. Split anchors each receive their own default window; one custom window cannot be applied to several split anchors.

```bash
gene-context stage1-select --results results/pou5_stage1 \
  --decisions /path/to/downloaded/loci.decisions.tsv \
  --out results/pou5_selected
```

Blank decisions, unknown IDs, duplicate rows, evidence changes and decisions from a different checkpoint block export. Outputs are `approved_loci.tsv`, `synteny_windows.bed`, the decision snapshot and provenance. BED has the standard six fields plus genome ID in column seven, since contig names are not unique across species. No context prediction starts automatically.

## Evidence and limits

The result contains the locus table, per-query support, every input genome's outcome, candidate DNA in source-locus orientation, complete raw BLAST and local miniprot GFF3 output, command receipts, source FASTA SHA-256 fingerprints and checksums of review evidence. FASTA hashing streams the input and adds an input-read pass. Raw BLAST files allow independent inspection. Work directories contain genome indexes/scratch data and are not themselves review evidence.

Manual `inspect` means investigate that interval. It does not prove orthology or functional pluripotency. Weak domains, paralogs, pseudogenes and assembly artifacts can survive. Lack of a match is not evidence of gene loss. The POU5 reference-window controls test mechanics; they do not validate sensitivity/specificity on all 18 complete assemblies or on unrelated gene families.
