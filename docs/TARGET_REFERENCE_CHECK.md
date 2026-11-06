# A small offline check before choosing target models

`check-targets` compares saved candidate proteins against a supplied local protein panel with BLASTP. It needs BLAST+ (`blastp` and `makeblastdb`), the candidate FASTA and the reference FASTA. Only the small protein panel is indexed, so BLASTP can use `--threads`; its direct `-subject` mode would ignore this setting. No genome indexing, gene prediction, model downloads, network access or GPU is involved.

Use this to support either review stop. At STOP 1, use an exported `evidence/<query_id>/protein.faa`. At STOP 2, check all coding options together:

```bash
gene-context check-targets \
  --proteins results/02_contexts/model_options.proteins.faa \
  --references references/target_and_paralogs.faa \
  --reference-labels references/families.tsv \
  --out results/02_reference_check --threads 4
```

Open `results/02_reference_check/report.html`, inspect the evidence, then record your choices in the original review decision files. This command never approves candidates or edits a checkpoint. Use a separate new output directory.

## Prepare the panel

Include trusted representatives of the target and related genes that could produce competing hits, preferably from several relevant lineages. A panel containing only the target cannot expose unrepresented alternatives. Gene symbols and existing database annotations also need scrutiny; this check does not resolve nomenclature.

The optional `families.tsv` must assign every reference FASTA ID exactly once:

```tsv
reference	family
chicken_POU5F3	POU5F3
zebrafish_POU5F3	POU5F3
human_POU3F2	POU3F2
```

Family labels are user-defined comparison groups. Here they group reference copies by gene name; they are not inferred evolutionary relationships. Without this file, reference names and coverage are reported but no family or competing-family score is inferred. `GN=` in a FASTA description supplies a display name only.

## Read the result

- `summary.tsv`: one row per candidate, its best reference, identity, coverage of the candidate **and** the reference, and any competing labelled family.
- `hits.tsv`: every reference with a significant match, ranked per candidate. Ties use reference ID for deterministic display; that ordering is not evidence favouring one tied hit.
- `blastp.tsv`: raw local alignments at E-value ≤ 1e-5, with no top-five truncation.
- `candidates.faa`, `references.faa`, `candidate_ids.tsv`, `reference_ids.tsv`: normalized search inputs and reversible IDs. Arbitrary input IDs are not passed directly to BLASTP.
- `manifest.json`: input/output SHA256 fingerprints, tool version, exact command and thresholds.

Coverage and identity are fractions from 0 to 1. Coverage counts aligned residues on each sequence separately; it is not alignment length divided by both lengths. Overlapping HSPs are greedily excluded on both axes before summing scores and weighted identities. This is a summary of local alignments, not a reconstructed global alignment.

Inspection flags:

| Flag | Meaning |
|---|---|
| `limited_coverage` | Best hit covers less than 50% of either protein. A short conserved domain can cover nearly all of a candidate fragment while covering little of a full reference. |
| `close_other_family` | Best hit to a different supplied family scores at least 95% of the overall best. Several references with the same family label do not count as competing families. |
| `internal_stop` | Candidate contains an internal stop; it was searched as X while keeping its position and this flag. |
| `no_reference_hit` | No significant match in this panel. This is not evidence of gene absence. |
| `reference_match` | At least one match, with none of the above flags. It is not automatic acceptance or proof of orthology. |

A single terminal stop is removed from FASTA inputs. References with internal stops and duplicate IDs are rejected. Gzip FASTA is accepted. DNA-only options have no protein to check and are absent from the input protein export.

The alternative score ratio is a comparison of local BLAST scores within this panel, not a probability. Sparse references, unequal evolutionary distances, missing paralogs and incomplete models can change the best hit. Use this evidence together with locus context, CDS structure and, when needed, a gene-family tree. No automatic orthology classifier is added by this command.
