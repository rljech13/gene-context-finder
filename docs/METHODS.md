# Methods and file reference

Start with the [README](../README.md) for the usual workflow. This page describes search rules, alignment settings and output fields.

## Genome manifest and annotations

`genomes.tsv` is tab-separated:

```tsv
genome_id	fasta	gff3
assembly_1	genomes/assembly_1.fna	annotations/assembly_1.gff3
assembly_2	genomes/assembly_2.fna.gz
```

Paths are relative to the manifest file, or absolute. `gff3` is optional. Use the annotation for the **exact assembly**, with matching FASTA contig IDs. FASTA and GFF3 can be gzipped. Convert GTF to GFF3 first. CDS Parent relationships define transcripts, including isoforms; gene features supply named neighbours. Missing annotation does not establish gene-free flanks.

Directory mode searches every FASTA recursively. Duplicate file stems require an explicit manifest with unique IDs. Directory mode does not guess annotation associations. Each file is an assembly observation: multiple assemblies of one species are not automatically merged into a copy count.

Inputs are read-only. Indices are created inside the run directory. Ordinary gzip inputs are decompressed there: allow space for uncompressed genomes and search indices. Use a new or empty output directory for each run. Stale results are not silently reused.

## What the pipeline does

1. **Protein (default):** TBLASTn searches the translated DNA of every contig independently of miniprot and GFF3. Defaults: E <=1e-5, >=30 aligned amino acids, >=30% identity. Individual exon seeds have no full-query coverage requirement. These local gates permit discovery through a conserved domain; they do not establish the gene family or orthology. Passing seeds on the same strand within `--dna-merge-distance` are collected as candidate regions. All passing regions remain in `dna_regions.tsv` with `source=tblastn`.
2. Overlapping same-strand annotated CDS models undergo an independent BLASTP sequence check using >=15% total query coverage (`--candidate-min-coverage`) as well as the length/identity/significance gates. A matching annotation supplies the full CDS even if the query matches only a domain. Its structural completeness and homology support are reported separately. GFF3 names never determine candidate acceptance.
3. For regions without passing annotated coding models, extract native-orientation DNA with `--refine-flank 50000` on either side and run miniprot locally with dense settings `-k3 -M0 -n1 -l3`. Only models whose coding exons overlap the original seed on the correct strand are considered; their translated proteins must also pass BLASTP. By default, at most `--max-refine-regions 20` regions per query/genome are refined in descending seed-score order. All unrefined or unresolved seeds remain visible; `protein_discovery/discovery.json` records the counts and cap. A resolved coding model may be partial or disrupted; sequence similarity does not repair missing exons automatically.
4. **Alternative fast protein route:** `--protein-search miniprot` maps proteins to the whole genome with the original >=50% query coverage and >=30% identity gates. Search/reporting controls include `--max-intron`, `--max-hits` and `--secondary-ratio`. Overlapping annotated proteins are checked with BLASTP. This route can miss divergent genes with short conserved domains, as demonstrated by the NANOG regression. New coding candidates from either protein route are rescored with BLASTP and ranked by summed non-overlapping HSP bitscore divided by the query's BLASTP self score. Raw discovery coverage/identity and the independent alignment coverage/identity are retained separately. `homology_support=local_alignment` explicitly marks a domain/fragment match; `extended_alignment` means the independent protein alignment also meets the full-coverage gate. Exact exon-coordinate duplicates are merged, preferring annotation; alternative models and separate loci remain visible.

**CDS:** translate the query and use the selected protein discovery route, and also run BLASTN. A terminal query stop is optional; internal stops and ambiguous coding triplets are rejected.

**DNA:** BLASTN searches all contigs on both strands. Passing HSPs overlapping annotated coding exons identify candidate transcripts. A complete unambiguous ATG-to-stop ORF in the query additionally enables the protein route, with that inference recorded. Otherwise, DNA without a matching coding annotation remains a DNA-region result. Arbitrary DNA does not uniquely determine CDS/protein. Protein-supported DNA seeds without a usable coding reconstruction also remain regions, labelled as translated-protein evidence rather than fabricated CDS models.

5. Re-extract CDS from source exons in transcript order and translate. Internal phases are checked without deleting bases from split codons. Separate terminal-stop features are handled explicitly. Missing start/stop gives a partial model; frameshifts, internal stops, ambiguous CDS and inconsistent phases give a disrupted model. All statuses remain in the inventory.
6. Extract gene+flank DNA and separate upstream/downstream windows in gene orientation. Flanks are anchored on the coding model, not a measured TSS. Record contig-boundary clipping. Retrieve annotated genes within the context radius and draw a native-coordinate context SVG.

DNA defaults: >=70% identity, E <=1e-5, HSP length >=30 bp. CDS associations need >=30 genomic CDS-overlap bases per HSP and >=10% union of query HSP spans (`--min-dna-coverage`). These are discovery criteria, not orthology proof. Same-orientation HSPs within 10 kb are merged only for region export; the region can contain several genes. Original BLAST/miniprot evidence is retained. Thresholds are operational and should be tuned to divergence, gene length, repeats and assembly quality.

## Alignments

MAFFT aligns proteins per query across qualifying candidates. CDS alignments are **backtranslated from protein alignments**, with exact codon/protein checks. Protein-only references have no invented reference CDS. Known query CDS supplies both. Terminal stops remain in CDS exports but are excluded from codon alignments.

Upstream/downstream windows are aligned separately in existing gene orientation. DNA uses progressive FFT-NS-2 by default (`--dna-alignment fast`); `--dna-alignment auto` allows MAFFT to choose potentially expensive iterative refinement. Protein alignment still uses MAFFT auto. Whole gene+flank alignment is opt-in (`--align-locus`) because introns and distant flanks may be long and poorly alignable. Complete predicted models enter alignments by default; `--include-partial` allows clean incomplete models. Disrupted models remain exported and are excluded. Arbitrary DNA queries also get DNA-region and region-flank alignments, even when no CDS model can be recovered. `--skip-alignments` runs search/extraction only.

Pairwise tables report identity and paired-position coverage. Gaps and ambiguous characters are excluded from comparable positions. Flank alignments are exploratory: identity and coverage must be considered together. An alignment does not demonstrate CRE activity, expression, individual-copy orthology or conserved function.

`--max-alignment-bases` bounds each input alignment set (default 2 million residues/bases). Exceeding it produces an explicit recorded error after extraction. Raise the limit or run extraction only; no arbitrary subset is silently selected.

## Repeat alignments without repeating genome searches

```bash
gene-context align --results results/cds_run --out results/cds_fast_alignment \
  --threads 8 --dna-alignment fast
```

This reads saved model/sequence exports, writes alignments to a separate directory,
and records the parent manifest and input hashes. It works after an alignment-stage
failure if extraction exports exist. The original result directory is unchanged;
source search errors remain recorded and are not cleared by successful realignment.

## Outputs

Open `report.html` first for a filterable overview of every input genome and links to recovered contexts.

| Path | Content |
| --- | --- |
| `report.html` | Filterable overview and context links |
| `genome_query_summary.tsv` | Every genome/query, including no hits and errors |
| `models.tsv` | Candidate transcripts, integrity, native coordinates, exon structure, evidence and flank completeness |
| `dna_regions.tsv` | DNA-hit envelopes, including matches without a coding model |
| `neighbors.tsv` | Annotated context genes for each coding model |
| `sequences/q0001/` | CDS, proteins, gene+flank DNA and separate oriented flanks across genomes |
| `genomes/<id>/models.gff3` | Coding models on native assembly coordinates |
| `genomes/<id>/*.context.svg` | Context overview; recovered model in blue |
| `genomes/<id>/dna_regions.fna`, `dna_regions.with_flanks.fna` | DNA-hit regions, with and without flanks, in query-relative orientation |
| `genomes/<id>/miniprot.gff3`, `blastn.tsv` | Original search evidence when relevant |
| `alignments/q0001/` | Protein, codon-CDS and flank alignments, pairwise comparisons |
| `manifest.json` | Input hashes, parameters, tool versions, commands and errors |

TSV coordinates are **0-based, half-open on the original genome**; GFF3 is 1-based, inclusive. FASTA headers link to model IDs. Zero-length flanks are omitted from FASTA and recorded as such in the model table.

A failed genome is recorded and the others are still processed. Exit codes: 1 for genome/alignment errors, 2 for invalid configuration, 0 for successful processing (including valid no-hit searches). No hit means not detected under these methods, not gene loss. Miniprot is less sensitive to distant homologs, and search reporting caps mean this is not a guaranteed exhaustive census.

## Query requirements

Supply one sequence type per FASTA file. Protein/CDS queries must encode at least 30 amino acids. CDS means spliced coding DNA in coding orientation, without UTRs; the terminal stop is optional. This release uses the standard nuclear genetic code. Query IDs are assigned in input order and original names are saved in `manifest.json`.

Each query is an independent search group. Cross-query hits are not merged into a gene family. Within each genome, models are ordered by query, contig and coordinates; this is not a ranking of orthology confidence.

## Tests and validation

Run `python -m unittest discover -s tests -v`. External-tool integration tests run when the required executables are installed and otherwise skip explicitly. See [validation](VALIDATION.md) and the [full-genome pilot](FULL_GENOME_PILOT.md).

Methods: [miniprot](https://github.com/lh3/miniprot), [NCBI BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/), [samtools faidx](https://www.htslib.org/doc/samtools-faidx.html), [MAFFT](https://mafft.cbrc.jp/alignment/software/).
