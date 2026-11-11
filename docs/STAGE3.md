# Stage 3: compare the reviewed models

This command accepts the **new locus-first stage2 checkpoint**. The older `compare-reviewed` command retains its legacy format.

## Review once, pass one JSON

In the stage-2 report, choose a target model, supply the region reason and save each region as proceed, excluded or unresolved. Only one representative model may be selected for a physical gene. Several distinct physical copies in one genome may share a comparison group and enter the same alignment. Copy identity is recorded separately in `copy_label` (native gene ID and coordinates by default); it does not assign orthology between species. Repeated observations of the same physical copy still cannot be selected twice. Use Add gene → Paralog / target copy for a supported additional copy in the same window. It uses its own model interval as the target seed; the other target-query, protein-family and coding-quality checks still apply. Other exceptional target choices require an explicit override and a comment. This records expert intent; it does not prove biological identity.

The main export downloads **one `stage2.decisions.json`**. Optional TSV tables are available separately. A portable draft is also accepted by the CLI only if all region decisions pass validation.

```bash
gene-context stage3 --results /mnt/project/nanog_stage2 \
  --decisions /path/to/stage2.decisions.json \
  --out /mnt/project/nanog_comparison --validate-only

gene-context stage3 --results /mnt/project/nanog_stage2 \
  --decisions /path/to/stage2.decisions.json \
  --out /mnt/project/nanog_comparison --threads 8
```

Validation checks the scientific seal, review ID, roles and region decisions, target-selection warnings/overrides, duplicate physical genes, comparison labels, neighboring-pair endpoints, CDS/protein hashes and backtranslation. When all CDS segments fit inside the saved genomic window, they are independently reconstructed from that DNA too. Edge models retain their limitations and require a review explanation; a selected partial target requires `--include-partial`. Disrupted CDS cannot enter comparison. Invalid decisions fail before creating an output directory or starting alignment tools.

## Computation

- Reuse the stored DNA, CDS and proteins. No genome download/search or Tiberius prediction is started.
- Run MAFFT protein alignment per explicit comparison label; add the supplied target panel from the saved stage-2 `target_support/subjects.faa`, where available. The BLAST query file contains extracted models and must not be used as the reference panel. Reference IDs are checked against `query_members` and recorded in the comparison manifest. References are excluded from the positional conservation profile. This profile counts selected gene copies; it is not a species-balanced conservation estimate.
- Backtranslate the protein alignment using the exact selected CDS, explicitly removing a terminal stop for alignment. Source CDS exports retain their terminal stop. Write pairwise protein/CDS identity, comparable length, coverage of both sequences and a positional protein conservation/coverage profile.
- Compare genomic windows as local BLASTN blocks (default), or optional LASTZ, on both strands. Write block coordinates, identity/coverage and coding/intronic/flanking summaries. Display blocks on the first selected locus in each comparison group; retain all pairwise tables.
- Compare only explicitly selected neighboring-gene models as final context evidence. Selected protein-pair proposals provide correspondence groups; components containing several physical copies in one species are not unique anchors. With no selected neighbors, the report states that context comparison has no reviewed anchors. It does not infer absence.

`--control-run` marks every result as a software control, used for validation with proposed test choices. It never upgrades them to user biological approvals.

## What to open and retrieve

Open `report.html`: a compact sequence reference sheet for sharing, with one row per selected copy, accession, native CDS/protein downloads, reference identity and both coverage directions. Click a comparison value to inspect and download its alignment. The earlier complete diagnostic view is retained as `technical.html`, with identity matrices, positional profiles, genomic blocks, full context and provenance.

Main exports:

- `sequence_catalogue.tsv`: target, species/assembly ID, copy label, physical locus, region/model IDs, original accession, source, native coordinates, coding quality, sequence hashes and recommendation reason.
- `selected.proteins.faa`, `selected.cds.fna`, `selected.genomic_windows.fna`: exact chosen inputs. Individual files under `sequences/` include genomic gene DNA only where the full model envelope fits inside the saved window; otherwise the catalogue records why it was not exported.
- `queries/group_*/proteins.aligned.faa`, `cds.codon.aligned.fna`, `protein_pairs.tsv`, `cds_pairs.tsv`, `protein_profile.tsv`: alignments and conservation measurements.
- `dna_pairs.tsv`, `dna_blocks.tsv`, `dna_categories.tsv`: genomic comparisons.
- `reviewed_neighbors.tsv`, `synteny.tsv`: explicitly chosen neighbors and their context summaries.
- `review_summary.tsv`: selected, excluded, unresolved and unsaved-stage-1 outcomes, preserving missingness.
- `decisions.snapshot.json`, `source.evidence.json`, `manifest.json`, `comparison.json`, `seal.json`: exact handoff and reproducibility evidence.

Identity always needs coverage and comparable length. Masking affects DNA search; genomic coverage denominators include all unambiguous A/C/G/T positions, including masked ones. A missing accepted block is not proof of biological loss. Comparison labels express the review scope: this command does not build a gene tree, establish copy orthology, verify expression/function, design an assay or identify experimentally active CREs. A taxon representative does not substitute for an unobserved organism's actual sequence.

On host2, run comparisons and tests through the existing constrained launcher, under the shared 250 GiB cap with zero swap, using `/mnt` for source snapshots, data, scratch and outputs. Keep source scientific checkpoints immutable.


## Compact reference sheet and technical appendix

The main sheet starts with the gene name and the selected-copy/reference counts, followed by one comparison table. Each supplied protein variant has its own column, including NANOG-like variants; no best-isoform choice silently replaces a query. Model details, alignment inspection, reference accessions and scope are collapsed by default. The report is entirely English, offline, responsive and printable.

Protein mode reuses the saved protein alignment and pair table against the supplied queries. CDS mode compares with selected models from those query species, where available. These can be different annotation variants, so the accession is explicit. Missing query CDS are not reconstructed from protein or downloaded implicitly. Self-references display 100% identity and coverage over comparable symbols, are labeled self, and are excluded from independent pair counts. An all-ambiguous sequence retains missing identity. Coverage is selected model / reference in both modes. The 80% coverage tint is a display cue only.

`sequence_reference_sheet.tsv` exports the per-copy sequence guidance and closest supplied protein reference (identity × smaller coverage). `reference_comparisons.tsv` exports all directional measurements and pair FASTA paths. Pair FASTAs are projections of the existing MSA with double-gap columns removed; the alignment is not optimized anew. The viewer renders the full displayed alignment with horizontal scrolling and always hides columns empty in every displayed row. Source MSA files retain all original columns; the original MSA span remains explicit. Full protein and codon-guided CDS MSA downloads remain available. Sequence guidance reports source/naming/completeness limits while retaining the exact user selection; it does not claim orthology, expression, function or interchangeability.

The technical appendix includes:


- A searchable sequence catalogue with original accessions, native copy labels, annotation/Tiberius CDS agreement, broad-database naming status, all retained query matches, and the exact chosen sequences.
- Protein and CDS identity matrices with coverage of both sequences and comparable length, a paged alignment viewer, and per-copy or species-balanced position profiles. The species-balanced profile gives each observed species one vote, divided equally among the observed copies at that alignment column.
- A separate within-genome copy table. Copy labels and short within-genome numbering do not assign orthology between species.
- Genomic alignment blocks against any selected reference copy, with native block coordinates and coding/intronic/upstream/downstream summaries. Overlapping windows from one assembly are flagged; shared physical DNA is not independent conservation evidence.
- A read-only context viewer that reuses the STOP-2 ruler, all available annotation/Tiberius coding alternatives, gene-level GFF context and saved protein correspondences. These descriptive proposals are explicitly separated from the manually selected anchors used by the final synteny statistics. Insufficient reviewed anchors remain unresolved.
- Command timings, source fingerprints, input decisions, checksums and downloadable tables/figures.

Final comparison categories and exported flanks use the selected coding-model envelope; the discovery seed is retained separately. These boundaries are not experimentally established transcript start/end sites. Upstream and downstream exports request 50 kb, stop at the saved window edges and record exact native coordinates, length and any clipping in `sequence_catalogue.tsv`.

New tables include `model_quality.tsv`, `target_pairs.tsv`, `paralog_pairs.tsv`, `conservation_summary.tsv` and `species_conservation_profile.tsv`. `selected.upstream.fna` and `selected.downstream.fna` and per-model flank files are exported. `START_HERE.html` and `README.md` make a complete result folder portable: extract it fully and open it in a modern browser with JavaScript; no internet or server is required for viewing.

Presentation can be refreshed without repeating scientific computation:

```bash
gene-context stage3-report --results COMPLETED_STAGE3 --out NEW_PRESENTATION
```

This verifies the completed-result seal, preserves alignments/sequences/decisions and regenerates the report from the stored data. Source rendering fingerprints are recorded separately in `presentation.json`.

## Final all-by-all comparison and model choices

The bottom of the main report contains the full protein/CDS alignment viewer, followed by the final model-choice table. Sequence A and B can be any stored sequence, including two native copies from the same genome or two protein query references. All sequences are shown together initially. An optional matrix provides navigation to each pair and its directional identity/coverage; self comparisons do not increase independent pair counts. Protein and CDS have separate input sets, so missing query CDS remain absent.

`all_alignment_pairs.tsv` indexes every stored independent pair and its FASTA under `all_alignments/`. Each pair is a projection of the saved MSA with double-gap columns removed, not a new pairwise optimization. Sequence identities, source alignments and decisions remain unchanged.

`primary_targets.tsv` (also `final_model_choices.tsv`) contains at most one proposed or reviewed target model per comparison gene/genome. `special_cases.tsv` contains the remaining STOP-2 selected copies and uncertain assignments. Both contain exact source identity, accession, assembly locus, selection reason, evidence and decision origin. Their CDS/protein FASTAs retain the exact native sequences. They partition only the selected STOP-2 set; previously discarded alternatives do not reappear.

The main page opens with all-by-all protein/CDS alignments, followed by these two tables. Change choice / Use as target opens a single editor, including a no-primary option. Saving requires a reason and marks that assignment reviewed. Existing primaries move to Separate cases when replaced. Notes for individual separate cases are editable. Changes persist per reviewed sequence set in browser storage. Download split package exports the **current** tables, four exact-sequence FASTAs and portable decisions JSON. Static files beside the report retain the generated state; browser edits do not write local files automatically.

Initial assignments are explicitly suggestions. Specifically known non-target copy labels are excluded from automatic primary selection. Remaining models rank by complete coding structure, then named target-query identity × minimum coverage; source identity and exact annotation/Tiberius CDS agreement break ties. If there is no query named for the target, the supplied panel is used explicitly. This heuristic does not establish orthology, gene function or family membership for unresolved copies. `gene_identities.tsv` retains source naming separately from comparison groups. Query preference is evidence, not automatic renaming.

```bash
gene-context stage3-finalize --results COMPLETED_STAGE3 \
  --decisions stage3.final.decisions.json --out NEW_FINAL_REPORT
```

This validates the exact source review and sequence-set fingerprint, every gene/genome bucket, permitted model IDs, reasons and reviewed flags. It generates a separate sealed report and split exports without running alignment or changing original STOP-2 selection. `finalization.decisions.json` and `finalization.evidence.json` retain assignment state and evidence. `model_choice_history.tsv` and `model_decision_audit.tsv` retain the original candidate and observation decisions in the technical appendix. No cross-species copy orthology is inferred from local copy numbering.
