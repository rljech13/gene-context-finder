# Workflow fixes and NANOG verification — 2026-11-09

Version 0.10.0. Open the [workflow entry page](http://127.0.0.1:8785/) or its local `index.html`. The working second review is under `stage2/`; `comparison_control/` is a separately marked software example. Original reports and user decisions were preserved. No Git commit or push was performed by this repair.

## What changed

| Original audit finding | Implemented change | Verification / remaining limit |
|---|---|---|
| A01: new STOP-2 has no comparison consumer | New `stage3` command accepts the sealed locus-first result and one decisions JSON; reuses saved models | Real MAFFT, codon backtranslation, BLASTN and catalogue exports; no discovery/prediction rerun |
| A02: a neighboring gene can silently become the target | Shared browser/backend target checks; explicit expert override plus reason for off-anchor/unsupported/conflicting-family choices | Negative controls reject the former SLC2A3-style mistake; documented expert intent remains possible |
| A03: window-dependent duplicate predictions | Predict once on each bounded union of overlapping windows; preserve original review windows and model observations | 25 windows mapped to 16 canonical intervals. 105 Tiberius observations reuse 70 unique predictions; 35 repeated predictions have zero CDS/protein conflicts |
| A03: duplicate or contradictory model choices | Reconcile physical genes conservatively; reject duplicate target selection and conflicting selected CDS for a physical gene | Distinct annotated genes and potential fusion models do not bridge separate physical loci. Separate copies require explicit comparison labels |
| A04: neighbor pairs survive excluded regions | Reject pairs whose region is excluded/unresolved or whose proposal IDs have changed | Shared UI/backend negative controls |
| A05: human NANOG is Unassigned | Fall back to the retained near-best symbols | Human prediction shows NANOG / NANOGP8; the remaining protein-level ambiguity is not hidden |
| A06: 18→13 scope is unclear | Show all input species and why they were not forwarded | Five first-review decisions were unsaved; this is not inferred gene absence |
| A07: actions/errors are hard to use | Quick stage-1 save controls, actionable species/coordinate errors, one primary STOP-2 JSON export and a concrete next command | Full 25-region isolated export/reload and real-browser presentation checks |
| A08: different windows look identical | Label displayed neighbor context separately from later prediction window | Original context and requested ±50 kb interpretation preserved |
| A09: sequence collection is separate | Per-model catalogue and exact protein/CDS/gene/window FASTA exports from reviewed choices | Assembly accession, source, native coordinates, hashes, quality and selection reason retained |
| A10: scale remains uncertain | Measure the real per-gene route and a dense 270-model UI control | The complete ~30-target ×18-species batch is not yet validated |

The emu second-copy observations now share the 284-aa model (previously 249/284); peafowl shares the 252-aa model (previously 252/271). This resolves inconsistent repeated runs, not biological correctness or orthology. All 11 proposed core-control CDS/proteins remained byte-identical to the earlier analyst proposal.

## Actual runs

- Stage 2: all 25 saved loci from 13 species; 16 CPU Tiberius intervals, 158 annotated model observations and 105 predicted observations; zero pipeline errors. Launcher wall time **1,229.31 s (20 min 29 s)**.
- Independent stage-2 audit: all seven checks passed. Original genomic FASTA/GFF hashes, all 25 window slices, all **263 native CDS/protein roundtrips**, raw prediction coordinates, GFF exports, broad BLASTP statistics and neighbor competition were checked. The native roundtrip check used BioPython independently of pipeline extraction.
- Corrected stage-3 software control: **11 proposed targets in 11 species**, 14 selected neighbor models, 11 selected pair proposals. **14.20 s** for comparison. Protein MSA: **16 records =11 chosen models +5 supplied references**, 659 columns; codon MSA: 11 records. Pairwise exports: 120 protein pairs including references, 55 CDS pairs, 55 genomic pairs and 1,624 accepted genomic blocks. Protein/CDS matrices display only chosen target models. Reference rows are excluded from the positional profile.
- Independent comparison audit passed: exact selected sequence hashes and translation, all protein/CDS pairwise numerators and denominators, source-reference identities, profile coverage, SVG parsing, catalogue files, all review outcomes and absence of search/prediction commands. All report links and seals were checked after local transfer.
- Final regression: 131 tests completed in 260.23 s: 130 passed, one optional LASTZ test skipped, zero failures. The final suite includes real cached CPU Tiberius and the adversarial target-reference integration check..
- Isolated DOM checks passed for all 25 loci, alternative sources, source-specific links, centering/scales, draft export/reload, wrong-checkpoint rejection and invalid choices. The dense 270-model fixture also passed. Test decisions stayed in isolated memory; working model choices remain pending.

A real-control run exposed an additional MSA input bug: extracted model proteins had been read from the BLAST query file as references. This was fixed to use the supplied panel in `target_support/subjects.faa`, checked against its query metadata. A regression fixture now deliberately supplies an unrelated protein in the BLAST query file. The initial 274-row control alignment is invalid and superseded; it is not linked from this package. The final independent audit requires the exact 16-record input set.

The initial independent-auditor launch lacked BioPython in system Python. It was rerun using the existing BioPython environment; no pipeline dependency was added and no scientific results were modified to make the audit pass.

## Using the result

1. Open `stage2/report.html`. Inspect Tiberius predictions and annotation alternatives, choose one exact target model per physical gene, choose useful neighbors/pairs, and save every region with a reason. Excluded or unresolved regions remain explicit.
2. Export one `stage2.decisions.json`. `stage3 --validate-only` checks the handoff without creating output or starting alignment tools.
3. Run `stage3` with that JSON and a new output directory. It writes alignments, pairwise conservation, reviewed-neighbor context and the exact sequence catalogue. Full commands and outputs: `docs/STAGE3.md` in the repository.

The first-review JSON is preserved under `stage2/selection.input.json`. Opening a different browser origin does not automatically import its old localStorage draft; use the report import control to reload it if revisiting STOP 1. The new second-stage evidence has its own review ID because predictions changed; old STOP-2 choices must not silently transfer.

## Remaining biological and batch work

The 11-model comparison is explicitly a software control, not a saved user biological approval, an ortholog set or a final recommended experimental sequence. The other reviewed regions remain unresolved in that control. Five species require completion of STOP 1 before a full 18-species STOP-2 run. Sequence/copy identity, uncertain CDS boundaries and suitable comparative representatives need the second manual review. Difficult copy relationships may need a focused family alignment/tree; no gene tree, expression test, experimental CRE claim or taxon proxy is inferred here.

A thin multi-target dispatcher/shared search-index cache and a measured ~30×18 run remain future work. Current checks cover NANOG and the existing test/control fixtures, not all target families or complete biological isoform recovery.

## Provenance and resource policy

Host2 jobs used the shared `gcfexperiments.slice`: aggregate hard cap **250 GiB (268,435,456,000 bytes)**, swap **0**. Source snapshots, inputs, caches, scratch and results are under `/mnt/data/i.rubinshteyn/gene_context_experiments`. Cumulative shared `oom` and `oom_kill` counters remained zero. Reported `max_child_rss_kib` is a process high-water mark, not aggregate job peak RAM.

Canonical producer snapshot archive SHA-256: `34dcfd35a602f39e28ae2d16c656e628ac56007bd96ab80e9af2d3a574dcd14b`. Final comparison/regression snapshot: `cb6848d4c15560582f015d6d95c366543ad879d22c394dd1a7d6e020743a7bcd`. The canonical prediction modules match the current checkout; the final runtime snapshot includes the corrected target-reference reader. Later documentation updates do not change runtime code.

Study-specific launch receipts, independent audits, control mapping, DOM checks, screenshots and helper-script hashes are in this output package’s `diagnostics/`. Scientific evidence is sealed separately from editable review decisions. The original audit remains a historical baseline, not a description of the repaired version.

## Initial display-window clarification (superseded below)

The stage-2 presentation now marks both physical display edges with native coordinates and the window length, and reports shown/total coding genes for each source. Isoforms share an arrow; noncoding annotations without CDS are not drawn. `//` marks models extending outside the displayed fragment; hidden groups are counted per edge. The collapsed details distinguish the displayed DNA from its shared canonical Tiberius input. In gene-order mode distances are schematic; in distance mode window widths use the shared bp scale.

For the chicken NANOG example, the displayed window is NC_052532.1:75826511–75931095 (104,584 bp), containing four annotation coding groups represented by six transcript models, versus three Tiberius predictions. The original coding annotations are SLC2A14, LOC101748883, NANOG and AICDA. This count does not describe all feature types in the complete source GFF.

This later presentation update changes the renderer, not the scientific outputs verified above. All scientific FASTA/GFF/evidence, decision TSVs and review IDs match the previous presentation. The 25-region DOM checks, negative decision controls and dense 270-model display checks pass; physical boundary orientation/counts and distance-scale proportions are included. See diagnostics/fragment_display_integrity.json and gcf_stage2_boundaries_* in the local output package.

## Combined annotation and Tiberius, including noncoding genes

This update supersedes the earlier coding-only presentation. The default diagram now combines both sources on separate arrows: annotation above, Tiberius below. Source-only genes remain visible. An ≠ badge means no identical CDS sequence is shared by an annotation/Tiberius alternative; it is not an automatic model rejection. Hollow annotation arrows represent gene-level GFF features without a CDS model in this window, including lncRNA. They have read-only context cards, not target-CDS decision controls or invented protein correspondences.

The actual original GFF files were reread for all 25 selected windows from 13 species; their SHA-256 values match the original stage-2 extraction receipts. The sealed `annotation_context.json` supplement contains 95 gene-feature observations (85 protein_coding, 10 lncRNA; overlapping windows repeat observations). Original evidence, all coding model FASTA/GFF files, window sequences, decision files and review IDs remain byte-identical.

Independent inspection of the original alligator GFF found five genes in the selected window: three protein-coding genes and lncRNAs LOC109281912 and LOC132250249. The chicken window contains six: four coding genes plus LOC121107507 and LOC112530152 (lncRNA). The alligator AICDA and FOXJ2 genes lie outside these display windows. The former three/four-gene presentation was omitting RNA annotations, not demonstrating their absence.

The shared annotation supplement was extracted on host2 in 12.67 s under the existing 250 GiB/zero-swap cap. Twenty-two focused Python checks passed on host2; the final 22-check local run also covers retaining the GFF supplement after protein reidentification. The combined-source/noncoding rendering, all 25 review windows, export/reload guards and the 270-model dense fixture passed DOM controls. The supplement is provided by `stage2-report --genomes ORIGINAL_MANIFEST`; new stage-2 runs retain gene features automatically.

The existing stage-3 software-control example remains its earlier explicitly chosen coding-neighbor comparison; these newly displayed RNA annotations do not become approved protein anchors. Current source/report files and measured receipts are separate from the original 131-test snapshot and canonical prediction run recorded above.
