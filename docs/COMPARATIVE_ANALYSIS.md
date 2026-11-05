# Comparative locus analysis

`analyze` accepts your sequences and your genome table. Gene names, species and reference assemblies are not built into the algorithm. Protein/CDS inputs use the standard nuclear genetic code; arbitrary DNA can also be compared as matched regions.

Since v0.3 the default context workflow predicts genes from the DNA of **every window**, then identifies the predicted proteins. See [window prediction](WINDOW_PREDICTION.md) for model selection, database preparation and outputs. The genome TSV must include `augustus_species`, or supply a global `--augustus-species MODEL`. The legacy v0.2 context workflow is explicitly selected with `--context-mode annotation`.

## One command

```bash
gene-context analyze --query targets.faa --query-type protein \
  --genomes genomes.tsv --out results/analysis --threads 8 --window 50000
```

Open `results/analysis/report.html`. Each target has a page with selected loci, gene-context diagrams, local DNA blocks, and protein/codon comparisons. The complete search inventory is retained under `search/`; comparative outputs are under `comparison/`.

The window includes the entire target plus the requested number of bases on each side. It uses the annotated gene envelope when the selected GFF3 model supplies one; otherwise it uses coding-model or DNA-hit boundaries. It is clipped at contig ends, oriented with the target, and records its native coordinates and SHA-256. Flanks do not establish a transcription start site or promoter boundary.

## Several references of one gene

By default, each FASTA record is a separate target group. To combine references, supply a tab-separated mapping covering every query exactly once:

```tsv
query	group
POU5F3_chicken	POU5F3
POU5F3_fish	POU5F3
NANOG_chicken	NANOG
```

Values in `query` must exactly match FASTA identifiers, or the recorded `q0001` IDs. Add `--query-groups groups.tsv` to the command. Matching loci are combined across those references; original queries and all candidate transcripts remain traceable.

## Candidate selection

Within each genome/target group, annotated isoforms of the same gene are grouped. Other models on the same strand sharing at least 80% of the smaller coding span are treated as alternative models of one candidate locus. Different annotated gene IDs remain distinct. The displayed representative is ranked by:

1. Query coverage multiplied by identity.
2. Complete coding model before partial or disrupted models when the sequence score ties; then coverage and protein length.
3. Model identifier for deterministic ties.

This is a sequence-based display choice, not a confirmed orthology assignment. `candidates.tsv` retains every alternative with its rank, query and selected flag. The overview explicitly reports multiple candidate loci. DNA-only fallback regions have no coding-model identity score and are ordered deterministically; inspect alternative regions before interpretation.

To override a display choice, use `--select-loci selection.tsv`:

```tsv
group	genome_id	model_id
POU5F3	my_assembly	my_assembly.q0001.m000003
```

Use exact values from the output; the group can be its label or internal ID. Omitted genome/group pairs retain automatic ranking. Invalid or duplicate selections are rejected.

## Neighbouring genes and local synteny

By default, AUGUSTUS predicts protein-coding gene structures from each window's DNA, with no protein or GFF3 hints. Predicted CDS and proteins are extracted, then BLASTP identifies the proteins against local/supplied references, or an automatically prepared Swiss-Prot library when no usable local reference proteins exist. Unnamed genes remain visible and can form homology groups across genomes. Provided GFF3 neighbours are saved separately for comparison and can supply naming reference proteins. Predictor model choice is explicit per genome or through `--augustus-species`.

With `--context-mode annotation` only, supplied GFF3 gene envelopes provide the primary neighbours. For windows lacking GFF3, BLASTX screens reference proteins, then miniprot reconstructs protein-supported coding models. This legacy route can use the same optional Swiss-Prot fallback, or `--context-reference none` for local inputs only.

Neighbour proteins are matched by BLASTP, independently of their annotation names. Default minimum identity is 30% and coverage is 50% of both proteins. Reciprocal unique best matches in the collected local/reference catalog form provisional homology groups. Their identity, coverage and reciprocal-match flags are exported. These are not genome-wide orthology calls.

Synteny compares gene order and strand **relative to the target**. At least two homology groups occurring exactly once in each window are required for order comparison. Repeated groups are excluded from the order score. The table gives the number of anchors, the fraction of concordant anchor pairs, and relative strand concordance. Too few anchors or missing annotation evidence is reported as insufficient/unknown. A narrow 50-kb window can contain too few genes even in a correct assembly; use a larger `--window` to investigate.

## Genomic DNA comparison

Every pair of selected windows in a target group is compared. The default engine is local BLASTN; `--dna-engine lastz` uses a separately installed LASTZ executable. LASTZ 1.04.60 was tested from the official source release. Choosing LASTZ makes it a required tool for that invocation.

Both engines search for local blocks rather than forcing the entire context into a global alignment. This preserves the distinction between aligned and unalignable DNA and allows reverse-oriented matches. Engine scoring differs: the same identity/block-length filters do not make the two engines equivalent sensitivity tests.

The pipeline preserves lowercase masking from the input genome. BLASTN uses lowercase seed masking and DUST; LASTZ uses its own masking/seeding rules. No repeat library is downloaded or inferred. Fully masked or highly divergent windows can yield no passing blocks; this is not zero identity or evidence of sequence loss.

After the engine-specific search, blocks need at least 50 paired canonical bases and 60% identity by default. These gates are controlled with `--dna-min-block` and `--dna-min-identity`. Higher-scoring blocks are accepted first; blocks overlapping an accepted block on either sequence are excluded from summary counts. This conservative filtering prevents duplicate coverage from repeated/overlapping alignments; raw blocks and rejection reasons remain available.

Identity counts exact matches among paired A/C/G/T bases, excluding gaps and ambiguity. Coverage is paired canonical bases divided by all canonical bases in that window. `dna_categories.tsv` separates coding, gaps between coding exons, target noncoding DNA and upstream/downstream DNA for each side. Between-CDS gaps are not a transcriptome-supported full intron/UTR annotation. Identity is blank when there is no passing alignment. Local matches, especially repeats, need review before claiming conserved regulatory elements.

The DNA plot uses `--reference-genome GENOME_ID`, or the first available selected locus. A requested reference with no available locus is reported as missing rather than silently replaced. Pairwise tables still cover all selected pairs. SVG files are vector graphics.

## Protein and CDS comparison

Complete selected models are aligned with MAFFT, together with the group's query proteins. Add `--include-partial` to include clean incomplete models as well; disrupted models remain excluded. Partial models can have the strongest sequence evidence and remain visible in the locus report even when excluded from the default protein alignment. CDS alignments are checked backtranslations of the protein alignment; a protein-only query contributes no invented reference CDS. Terminal stops are retained in CDS exports and omitted from codon alignments.

The page shows pairwise amino-acid identity/coverage and a position profile. Query references are excluded from the profile, so adding many query proteins does not increase apparent cross-species support. Consensus fractions require at least two observed locus residues; the profile also reports locus coverage. These descriptive comparisons are not selection tests or phylogenetic conservation scores.

## Reuse saved searches

```bash
gene-context compare --results results/previous_search \
  --out results/new_comparison --threads 8 --window 50000 \
  --augustus-species YOUR_MODEL
```

Original searches and their manifests remain unchanged. The source genomes (or saved working FASTA caches) and supplied GFF3 files must remain accessible on the machine running comparison. Target sequence roundtrips and GFF3 hashes are checked, and each newly extracted window is fingerprinted. Copying only the original HTML report is not enough to recompute context.

The provenance manifest records original search errors as well as comparison errors. A successful comparison does not erase a source-search failure. Missing biological evidence is distinct from a failed executable/input. Use separate new/empty output directories for each run. The default maximum window is two million bases; raise `--max-window-bases` explicitly for longer genes. Protein alignment input is bounded by `--max-alignment-bases`.
