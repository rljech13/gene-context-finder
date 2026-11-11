"""Portable, offline report for reviewed sequence conservation and copy identity."""
from collections import Counter
import json
from pathlib import Path

from .comparison_report import dna_svg, _protein_svg
from .sequence import read_fasta
from .stage3_analysis import read_table


def refresh(args):
    import shutil
    from .sequence import sha256
    from .stage1 import new_output,dump
    from .stage3_context import write_context_page
    source=Path(args.results).resolve();seal=json.loads((source/'seal.json').read_text())
    for name,digest in seal.items():
        p=(source/name).resolve()
        if source not in p.parents or not p.is_file() or sha256(p)!=digest:raise ValueError('Stage-3 result changed: '+name)
    manifest=json.loads((source/'manifest.json').read_text())
    if manifest.get('stage')!='stage3_comparison':raise ValueError('Not a stage-3 result')
    final_state=None
    if getattr(args,'command','')=='stage3-finalize':
        from .stage3_finalization import validate_choices
        evidence=source/'finalization.evidence.json'
        if not evidence.exists():raise ValueError('Refresh this result with stage3-report before final assignment')
        final_state=validate_choices(json.loads(Path(args.decisions).read_text()),json.loads(evidence.read_text())['spec'])
    target=Path(args.out).resolve()
    if target==source or target in source.parents or source in target.parents:raise ValueError('Use a separate presentation directory')
    root=new_output(target)
    for name in seal:
        p=root/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source/name,p)
    if final_state is not None:dump(root/'finalization.decisions.json',final_state)
    comparison=json.loads((root/'comparison.json').read_text());data=json.loads((root/'source.evidence.json').read_text())
    if (root/'context_evidence.json').exists():write_context_page(root,json.loads((root/'context_evidence.json').read_text()))
    write_report(root,data,read_table(root/'sequence_catalogue.tsv'),comparison['groups'],comparison['selected'],read_table(root/'review_summary.tsv'),comparison.get('blocks',[]),comparison['context'],comparison['errors'],manifest['control_run'],comparison['analysis'],manifest)
    dump(root/'presentation.json',dict(source=str(source),source_manifest_sha256=sha256(source/'manifest.json'),note='Stored scientific alignments, sequences and STOP-2 decisions retained. Final target/special-case assignment is stored separately; suggestions are not manual approvals.',finalization_decisions_sha256=sha256(root/'finalization.decisions.json'),renderer_sha256={name:sha256(Path(__file__).with_name(name)) for name in ['stage3_report.py','stage3_report.txt','stage3_summary.py','stage3_summary.txt','stage3_handoff.py','stage3_finalization.py','stage3_finalization.txt','gene_identity.py','stage3_context.py','comparison_report.py','stage2_neighborhood.txt','stage2_cluster.txt']}))
    dump(root/'seal.json',{str(p.relative_to(root)):sha256(p) for p in root.rglob('*') if p.is_file() and p.name!='seal.json'})
    print('Updated stage-3 report: '+str(root/'report.html'));return 0


def write_report(root,data,catalogue,groups,selected,observations,blocks,context,errors,control,analysis,manifest):
    root=Path(root)
    catalogue=[dict(row,**{k:int(row[k]) for k in ['protein_aa','start','end','upstream_bp','downstream_bp'] if row.get(k) not in ('',None)}) for row in catalogue]
    observations=[dict(row,compared=row['compared'] if isinstance(row['compared'],bool) else str(row['compared']).lower()=='true') for row in observations]
    quality={r['model_id']:r for r in analysis['quality']};models=[];alignments={}
    counts=Counter((r['group'],r['genome_id']) for r in selected);indices=Counter()
    for row in selected:
        key=row['group'],row['genome_id'];indices[key]+=1
        species=row['genome_id'].split('_GC')[0].replace('_',' ');words=species.split();short=(words[0][0]+'. '+words[1]) if len(words)>1 else species
        short+=(' · '+str(indices[key])) if counts[key]>1 else ''
        item=next(c for c in catalogue if c['model_id']==row['model_id'])
        models.append(dict(item,group=row['group'],short_label=short,quality=quality[row['model_id']],target_support=row['target_support'],identification=row['identification']))
    label_by_id={m['model_id']:m['short_label'] for m in models}
    for group in groups:
        folder=root/'queries'/group['id'];chosen=[dict(r,display_label=label_by_id[r['model_id']]) for r in selected if r['group']==group['id']]
        profiles=[dict(coverage=r['species_coverage'],consensus_fraction=r['species_balanced_consensus_fraction']) for r in analysis['species_profiles'] if r['group']==group['id']]
        if profiles:_protein_svg(folder/'species_conservation.svg',profiles,title='Species-balanced protein conservation (query references excluded)',unit='species')
        for row in chosen:dna_svg(folder/('genomic_blocks.'+row['model_id']+'.svg'),group['name'],chosen,blocks,row['model_id'])
        alignments[group['id']]={kind:read_fasta(folder/name) if (folder/name).exists() else {} for kind,name in [('protein','proteins.aligned.faa'),('cds','cds.codon.aligned.fna')]}
    payload=dict(review_id=data['review_id'],groups=groups,models=models,observations=observations,
        summary=analysis['groups'],pairs=analysis['pairs'],paralogs=analysis['paralogs'],dna=read_table(root/'dna_pairs.tsv'),
        categories=read_table(root/'dna_categories.tsv'),synteny=read_table(root/'synteny.tsv'),
        context_models=len(context),reviewed_pair_count=analysis['reviewed_cross_locus_pairs'],
        alignments=alignments,query_members=data.get('query_members',{}),errors=errors,control=control,
        naming_status=analysis['naming_status'],source_agreement=analysis['source_agreement'],
        overlapping_window_pairs=analysis['overlapping_window_pairs'],parameters=manifest.get('parameters',{}),
        timings=manifest.get('timings',{}),core_elapsed_seconds=manifest.get('elapsed_seconds',0),resources=manifest.get('resources',{}))
    text=json.dumps(payload,ensure_ascii=True,separators=(',',':')).replace('<','\\u003c').replace('&','\\u0026')
    template=Path(__file__).with_name('stage3_report.txt').read_text()
    (root/'technical.html').write_text(template.replace('__DATA__',text))
    from .stage3_summary import write_summary
    write_summary(root,payload)
    (root/'START_HERE.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Gene sequence reference sheet</title><style>body{font:17px/1.6 system-ui;max-width:760px;margin:50px auto;padding:24px;color:#213a3d;background:#f3f6f5}a{color:#1d686b}h1{line-height:1.2}</style><h1>Gene sequence reference sheet</h1><p><a href="report.html">Open the reference sheet →</a></p><p>The report opens with complete protein / CDS alignments and any-pair comparisons. Two final tables separate one proposed target representative per organism from additional copies and uncertain cases. Each row explains the assignment; you can change it and export the split sequence set. Reference-comparison tables and diagnostics are in the technical appendix.</p><p>Extract the entire archive before opening this page. Use a modern browser with JavaScript enabled. The report works offline; keep its folders together. No installation or server is needed.</p><p><a href="technical.html">Technical appendix</a> · <a href="README.md">Instructions</a> · <a href="sequence_reference_sheet.tsv">Sequence choices table</a></p></html>')
    (root/'README.md').write_text('''# Reviewed sequence conservation report

Open START_HERE.html, then the reference sheet in a modern browser with JavaScript enabled. Keep the entire folder together. No server, Python or internet connection is needed to view these results.

report.html opens with the complete all-by-all protein / CDS alignment viewer and ends with two editable assignment tables: one target representative per gene/genome and separate cases, with native sequence downloads. Select any two sequences, or use the matrix, to inspect and download a pair. The former initial reference-comparison table is omitted from the main page; its measurements remain in reference_comparisons.tsv and the technical appendix. The entire displayed alignment is rendered with horizontal scrolling; columns empty in every displayed row are always hidden. Source MSA downloads retain the original columns. The final all-by-all block includes all selected copies and supplied protein references, any ordered sequence pair, an identity/coverage matrix, and every pair FASTA. The final assignment follows this block and partitions the exact STOP-2 selected native models into Target gene and Separate cases, with accessions, genomic loci, assignment reasons and sequence downloads. The all-model protein and CDS alignments remain directly downloadable. Technical figures, context and provenance are in technical.html and the supporting files.

Protein alignments include selected models and the supplied query proteins; CDS alignments include the selected models whose CDS are available. Sequence selectors identify each accession. No missing query CDS is invented or retrieved. A reference compared with itself displays 100% identity and coverage over comparable symbols, with a self label; it is excluded from independent pair counts. With no comparable symbols, identity remains missing. Missing measurements remain explicit. Coverage is always selected sequence / reference. Always consider identity together with the coverage of both sequences.

sequence_reference_sheet.tsv contains the sequence guidance and closest supplied protein reference (ranked by identity times the smaller coverage). reference_comparisons.tsv contains the exact comparable positions, identity, both coverage directions and pair-alignment paths. Pair downloads are projections of the stored multiple alignment with double-gap columns removed, not new pairwise optimization. A close match does not establish orthology or functional interchangeability. Sequence guidance retains the exact manual choice and reports integrity, source disagreement and naming limits; it is not experimental validation.

The report includes selected source sequences and model quality, protein and codon-guided CDS alignments, identity and coverage matrices, within-genome copy comparisons, genomic local-alignment blocks, exact upstream/downstream exports, and the available annotation/Tiberius gene context. All tables, alignments, figures and provenance are linked.

Related copies may share a comparison group. Copy numbers in short display labels are local to each genome and do not assign cross-species orthology. Query references guide protein alignment and are excluded from the selected-copy conservation profile. The species-balanced profile gives each observed species one vote, divided equally among its observed copies at that column.

Genomic DNA is compared as local blocks, not a forced global alignment. DNA matrices require identity, coverage and comparable length together. Overlapping windows from the same assembly are flagged: similarity may reuse the same physical DNA, rather than independently support conserved paralogous regions. Model-span DNA and flanks use the selected coding-model envelope, not an experimentally established transcription start/end site. Flanks stop at the saved window edges.

The context overview retains all available models and saved protein-match proposals from STOP 2, including unreviewed and noncoding context. Final synteny statistics use only manually selected neighboring models and correspondence pairs. Too few reviewed anchors means synteny is unresolved; it does not mean genes are absent. A conserved block does not establish CRE activity, and sequence comparison does not establish expression or experimental function.

The exact selected models are user choices. Quality flags identify unresolved naming, disagreements between annotation and prediction, and partial/edge or short-support cases; they do not silently substitute a new model. Inspect these flags before choosing a sequence for experimental work.

all_alignment_pairs.tsv indexes all independent protein/CDS pairs, with both coverage directions and paths under all_alignments/. These files are MSA projections, not newly optimized pairwise alignments. final_model_choices.tsv is an alias of primary_targets.tsv (one proposed or reviewed target representative per gene/genome). special_cases.tsv retains additional copies and uncertain assignments. model_choice_history.tsv retains the earlier candidate outcomes separately in the technical appendix. model_decision_audit.tsv retains every source model observation, exact saved model/region decisions and comments. Context-only neighbors are separate from target candidates. An unselected alternative is not a biological rejection; unresolved cases remain unresolved.

Final assignment is separate from the STOP-2 model selection. Initial suggestions rank complete coding models first, then target-query protein identity times the smaller model/reference coverage. Exact source naming and annotation/Tiberius CDS agreement break ties. Known specifically named non-target copies are never automatically renamed as the target. With no matching named target query, the supplied reference panel is used explicitly. This ranking is a proposal, not resolved orthology or functional validation. Exact copy identity is exported separately in gene_identities.tsv; a family-level BLAST assignment does not establish a specific paralog name.

Use Change choice or Use as target to choose one representative, or no representative, for an organism. Enter the reason and click Save reviewed choice. The former primary is then a separate case. Edit note records a separate-case explanation. Changes stay in this browser under this exact reviewed sequence set; Download split package creates current primary/special-case TSV tables, exact native CDS/protein FASTAs and stage3.final.decisions.json. Export/import choices transfers this state between browsers. Static TSV/FASTA files alongside the report reflect the state when the report was generated, so use the ZIP button after browser edits. Proposals and manually reviewed assignments are identified separately. TSV coordinates are zero-based inclusive start, exclusive end; displayed genomic positions are one-based inclusive.

To materialize exported choices as another sealed full report without recomputing alignments:

gene-context stage3-finalize --results COMPLETED_STAGE3 --decisions stage3.final.decisions.json --out NEW_FINAL_REPORT

The command checks the exact review and selected sequence set, one valid model per gene/genome, and a reason for each assignment. Original STOP-2 decisions, exact sequences and all-by-all alignments are retained. Selected copies moved to Separate cases are not biologically rejected and remain in the all-by-all alignments. Neither table reintroduces models discarded at STOP 2. finalization.evidence.json records selection evidence and grouping; finalization.decisions.json stores the presentation's assignment state.

Reproducibility: decisions.snapshot.json, source.evidence.json, manifest.json, comparison.json and seal.json retain inputs, decisions, settings, command timings and file checksums.
''')
