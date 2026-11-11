# Work plan: target-gene conservation and sequence selection

Updated for version 0.10.0 on 2026-11-09. The [original audit](WORKFLOW_AUDIT_2026-11-09.md) is a historical baseline: it records the reproduced gaps before these fixes. The current implementation includes shared prediction intervals, shared STOP-2 validation, and the new `stage3` comparison consumer. Study-specific run receipts distinguish software controls from biological review.

## Intended result

For each requested protein-coding gene, describe conservation across the supplied species at three levels: protein, coding DNA, and surrounding genomic DNA. Provide a traceable, recommended sequence/model for work in each organism, together with alternatives and an explicit confidence/reason. Taxon-level recommendations identify a representative species and the supporting species; they are not substituted for an unobserved organism's actual sequence.

The present study is approximately 30 targets by 18 species (540 target/species cells, potentially more physical loci and transcript models). The core analysis follows the requested gene. Additional family members are recorded as possible paralogs and retained for later expansion. A family-wide census, functional equivalence and expression validation are outside this first deliverable. A target with unresolved copy identity remains unresolved rather than being forced into the core set.

## Keep the existing execution structure

1. **Stage 1 / STOP 1:** multi-reference search, candidate-locus ranking, per-query evidence and annotation-based preliminary context. Review and select physical loci/windows.
2. **Stage 2 / STOP 2:** extract DNA; retain annotated alternatives; predict coding genes with Tiberius; identify predicted proteins; compare neighbors; review the exact target CDS/protein models and proposed correspondences. Annotation-only models may be chosen explicitly when supported; prediction is independent evidence rather than a required replacement of an annotation.
3. **Final comparison:** consume the saved STOP-2 selections, compare the chosen corresponding loci/models, and generate the conservation report and recommended-sequence catalogue. Do not rerun discovery or prediction merely to perform comparison.

Retain the two existing manual review gates. Display controls and diagnostic suggestions do not constitute approvals. Ambiguous cases can stop without blocking other target/species cells.

## Existing implementation and actual remaining work

| Component | Current coverage | Required work |
|---|---|---|
| Multi-reference discovery and ranked candidate loci | Implemented in `stage1.py`, `stage1_proteins.py`, `locus_chaining.py`, `locus_ranking.py`; BLASTP over annotation products, explicit genome rescue, per-query support, editable filters | Preserve this route; verify reference-panel metadata and retain original symbols/isoforms. Automatic rescue triggers should be explicit, not assumed to exist. |
| Preliminary annotated synteny and STOP 1 | Implemented in `early_synteny.py`, stage-1 context/report modules | Polish interpretation and reuse source context; keep target similarity separate from independent neighbor evidence. |
| DNA extraction, Tiberius, broad protein identification | Implemented in `stage2.py`, `tiberius_backend.py`, `stage2_db.py`; complete NANOG run reached STOP 2 for 25 saved loci across 13 species | Shared bounded prediction intervals are implemented in `stage2_intervals.py`; original review windows retain their mappings. Predictions are reused across observations. Verify exact model agreement and scientific output on the canonical rerun. This is not validation of all 18 species or 30 targets. |
| Stage-2 review and exports | Implemented in `stage2_report.py`, `stage2_review.txt`, `stage2_neighborhood.txt`; Tiberius-only diagram is the default, annotation/combined views remain available | Shared backend/UI constraints reject unsupported target selection without documented override, duplicate physical targets, inconsistent selected CDS and contradictory region/pair choices. Scope, source, alternative names and actionable errors are now visible. Preserve scientific evidence and decisions when refreshing presentation. |
| Protein MSA, codon-guided CDS alignment and pairwise statistics | Existing reusable code in `alignment.py` and `comparison_report.py`; legacy route also writes a positional protein profile | Implemented in `stage3.py`: validated sealed STOP-2 JSON handoff, explicit comparison groups, exact model provenance, native-DNA/backtranslation checks, terminal-stop handling, pairwise identity/coverage and positional profile. Domain-specific interpretation remains future work. |
| Local genomic DNA comparisons and context summary | Existing legacy code in `dna_compare.py`, `context_compare.py`, `comparison.py`: local blocks, coverage, identity, coding/intron/flank categories and neighbor order | Reused by `stage3`: local blocks, coverage, coding/intron/flank categories, explicitly selected neighbor context. Masking and clipping limitations are stated. No forced whole-window MSA. |
| Recommendation of sequences for an organism/taxon | Model alternatives and provenance exist; an analyst proposal was exported for NANOG | Implemented per-model catalogue and protein/CDS/gene/window FASTA exports tied to exact review choices, hashes, assembly and reasons. Biological acceptance and taxon proxy justification remain the reviewer’s task; no automatic longest-model rule. |
| Multi-target batch orchestration | Per-target runs and stage-1 per-genome workers exist; annotation catalogs/reference weights can be shared | Add a thin target manifest, checkpoint-aware dispatcher and one project summary. Reuse safe caches and bound aggregate resource use. Stage-1 BLAST databases are currently rebuilt in each search work directory; shared index reuse remains work, not an existing optimization. |

`compare-reviewed` retains its legacy model-first format. The new bridge is the separate [`stage3` command](STAGE3.md), consuming the locus-first stage-2 checkpoint and its one authoritative decisions JSON. It does not repeat discovery or prediction.

## Priority 1: stable loci and consistent models

- Create one canonical prediction interval for overlapping selected neighborhoods on a given assembly/contig. Preserve each selected locus and its original requested display window as a mapping onto that interval; avoid unbounded merging of a chain of overlapping windows. Respect the existing maximum window/resource limits and flag intervals that cannot be processed consistently within them.
- Predict once per canonical interval/configuration and reuse the result. Cache keys include the assembly/input hashes, interval, strand convention, masking, predictor/model versions and parameters. This removes duplicate context-dependent runs; it does not guarantee that a model is biologically correct.
- Separate **physical locus**, **coding-model alternative**, **repeated observation**, and **possible paralog**. Use assembly coordinates, strand and compatible exon structure as evidence for reconciliation; mere overlap cannot safely merge distinct genes. Ambiguous fusion/split assignments stay flagged.
- Compare annotation and prediction at exon/CDS level. Record exact agreement, altered starts/stops, splice differences, partial/edge models and missing prediction separately. Resolve the emu 249/284-aa and peafowl 252/271-aa examples before selecting their final models. Existing prediction isolation tests cover unrelated input-window effects, not invariance to shifting the same window boundaries.

Completion criterion: the same physical gene reached through overlapping selected windows is represented once with traceable model alternatives; every original selected window maps back to the common prediction evidence.

## Priority 2: select the core target and sequence

Use the requested reference identity, locus/annotation provenance, coding integrity, protein support and context together. Store the chosen model ID, sequence source, CDS/protein hashes, exact assembly coordinates and a short reason. Retain unchosen annotated and predicted alternatives.

The core comparison set can include a supported annotated model even when Tiberius misses or truncates it; the source must be explicit. Structural completeness alone establishes neither expression nor protein function. Unsupported models do not become complete by selecting a longer alternative.

Check competing local copies before accepting a one-to-one comparison. For an unambiguous supported target, proceed after review. When copy relationships remain uncertain, use a focused diagnostic family alignment/tree plus the species relationships and context, or defer that target/species cell with `copy_unresolved`. This targeted check is not a requirement to characterize every member of every family now.

Completion criterion: each target/species cell has either a supported selected model, an explicit retained set of alternatives with its comparison scope, or a reason it is not ready. The decision is reproducible from saved evidence.

## Priority 3: connect the existing comparison modules

Read the sealed stage-2 evidence and the exported decisions. Validate IDs, roles, assembly/model/sequence hashes and contradictions. Reuse stored DNA/CDS/proteins. Never replace the selected model silently or start prediction again.

For each approved comparison group:

- **Protein:** MSA, pairwise identity and coverage in both directions, comparable-position counts, a positional conservation/coverage profile, gaps and uncertain/poorly alignable regions. Show whole-protein and supported domain/region summaries where annotations are available.
- **CDS:** codon-guided alignment from the same selected protein models, checked backtranslation and nucleotide identity/coverage. Preserve terminal stops in source exports and handle them explicitly in the alignment. Do not add dN/dS or selection claims to this first release.
- **Genomic context:** local homologous blocks with native coordinates, identity and covered fractions of each region, separated into coding, intronic and flanking context. Report no detected alignable blocks and technical missingness separately. Conservation does not establish CRE activity.

Always report identity together with comparable length and coverage. A small conserved domain must not look like a conserved full gene. Missing or unreliable positions are not automatically divergent positions. The report must support taxon-specific comparisons rather than averaging all 18 species into a universal conserved/nonconserved label. Gene rankings, if offered, require comparable species sets and coverage criteria.

Completion criterion: a small saved STOP-2 selection reaches all three comparison outputs without discovery/prediction reruns, and each aligned record traces to the exact approved model.

## Priority 4: report and recommended-sequence catalogue

Default overview: one row per target/species cell, selected copy/model and source, coding quality, context support, comparison readiness and recommendation reason. Each gene page shows compact protein/CDS comparisons, position profiles, genomic blocks and links to full evidence. Retain no-hit, failure, partial and unresolved cells.

Export a sequence catalogue with target, species/taxon, assembly accession, locus/model/transcript/protein accessions and versions, source (`annotation`/`Tiberius`), coordinates, hashes, quality flags, review status and recommendation reason. Export the matching protein, CDS and genomic sequences. Keep biological isoforms as labelled alternatives; a representative for comparative analysis is not automatically the right transcript for every experiment. Name the organism and limitations of any taxon-level proxy. No primer, construct or validated regulatory-element design is implied.

Completion criterion: the user can see what is conserved, which regions are uncertain/divergent, and retrieve the exact recommended sequence for a chosen organism without reconstructing it from intermediate files.

## Priority 5: scale the proven path

Use a target manifest to dispatch the same existing stages for approximately 30 panels and 18 fixed assemblies. Keep target runs separate while sharing content-verified catalogs, indexes, reference databases and model weights where supported. Resume from checkpoints and regenerate reports without recomputing scientific evidence. Shared-cache creation requires collision protection and fingerprint checks.

Keep host2 experiments under the existing aggregate 250 GiB RAM cap and zero swap, with inputs, temporary data and outputs under `/mnt`. Record wall time, CPU/RAM measurements, tool versions and scientific error states. Failures in one target/species cell must remain visible without discarding valid results elsewhere.

Validation is incremental: use the current NANOG case to close the known locus/model issues, then a small set covering a simpler target, close competing copies and a difficult coding structure. Verify the new handoff, exact sequence roundtrips, comparison denominators and manual-decision preservation. Finally run the complete 18-species panel and scale to the remaining targets. Existing audits are reusable evidence; they do not substitute for testing the new adapter or the complete batch.

## Deferred work

Complete gene-family evolution, exhaustive paralog discovery, RNA-supported isoform validation, expression-stage/function claims, selection tests, experimentally validated CREs and assay/construct design are separate extensions. Preserve the evidence needed to revisit them without making them prerequisites for every core target.
