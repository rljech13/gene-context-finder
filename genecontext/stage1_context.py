"""Attach automatic annotation-based context evidence to locus-first reviews."""
from pathlib import Path
import shutil
import sys
from uuid import uuid4

from .sequence import write_tsv, sha256
from .tools import Runner
from .stage1 import dump, new_output, now, verify, write_outputs


def attach(args, root, data, runner=None):
    if getattr(args,'early_synteny','auto')=='off':
        return False
    from .early_synteny import run_early_synteny, settings
    runner=runner or Runner(root)
    original=list(data.get('commands',[]));start=len(runner.commands)
    manifest=getattr(args,'datasets',None) or getattr(args,'genomes',None)
    if not manifest:
        manifest=root/'datasets.input.tsv'
        fields=['genome_id','fasta','gff3','protein','cds','annotation_provenance','species','assembly_accession']
        write_tsv(manifest,data['genomes'],fields)
    cache=getattr(args,'neighbor_cache',None) or root.parent/'_genecontext_cache'/'early_synteny'
    output=root/'raw'/'early_synteny'
    print('[stage 1] automatic early synteny',file=sys.stderr,flush=True)
    try:
        inputs={}
        receipts=root/'raw'/('context_inputs_'+data['review_id'])
        receipts.mkdir()
        for name,value in [('datasets',manifest),('source_datasets',getattr(args,'source_datasets',None))]:
            if value:
                path=Path(value).resolve()
                digest=sha256(path)
                shutil.copyfile(path,receipts/(name+'.tsv'))
                if sha256(receipts/(name+'.tsv'))!=digest:raise ValueError('Dataset manifest changed while copying: '+str(path))
                inputs[name]=dict(path=str(path),sha256=digest,snapshot=str((receipts/(name+'.tsv')).relative_to(root)))
        result=run_early_synteny(data,genomes_manifest=manifest,
            source_genomes_manifest=getattr(args,'source_datasets',None),cache_dir=cache,
            out_dir=output,threads=args.threads,config=settings(args),runner=runner)
        known={(l['genome_id'],l['locus_id']) for l in data['loci']};seen=set()
        for item in result['loci']:
            key=(item['genome_id'],item['locus_id'])
            if key not in known or key in seen:raise ValueError('Unknown/duplicate locus in early-synteny evidence: '+str(key))
            seen.add(key)
        if seen!=known:raise ValueError('Early-synteny evidence has missing loci')
        for item in inputs.values():
            if sha256(item['path'])!=item['sha256']:raise ValueError('Dataset manifest changed during context analysis: '+item['path'])
        result['inputs']=inputs
    except (OSError,ValueError,RuntimeError) as error:
        print('[early synteny] '+str(error),file=sys.stderr)
        result=dict(schema='early_synteny_v1',status='failed',parameters=settings(args),references=[],source_resolution=[],datasets=[],
                    errors=[dict(stage='early_synteny',error=str(error))],loci=[dict(genome_id=l['genome_id'],locus_id=l['locus_id'],
                    status='early_synteny_error',neighbors_requested=getattr(args,'synteny_neighbors',5),neighbors_used=0,
                    window_start=l['start'],window_end=l['end'],neighbors=[],comparisons=[],flags=[str(error)]) for l in data['loci']])
    result['artifact']='raw/early_synteny/evidence.json'
    result['source_review_id']=data.get('parent_review',data['review_id'])
    output.mkdir(parents=True,exist_ok=True);dump(output/'evidence.json',result)
    data['early_synteny']={k:v for k,v in result.items() if k!='loci'}
    mapping={(l['genome_id'],l['locus_id']):l for l in result['loci']}
    for locus in data['loci']:locus['early_synteny']=mapping[locus['genome_id'],locus['locus_id']]
    data['commands']=original+runner.commands[start:]
    return bool(result.get('errors'))


def context(args):
    source,data=verify(args.results)
    root=new_output(args.out)
    if source in root.parents or root in source.parents:raise ValueError('Context report must be separate from its source')
    for name in ['queries.fasta','loci.fna']:
        if (source/name).exists():shutil.copyfile(source/name,root/name)
    shutil.copytree(source/'raw',root/'raw',ignore=lambda folder,names:['early_synteny'] if Path(folder)==source/'raw' and 'early_synteny' in names else [])
    if (source/'raw/early_synteny').exists():
        shutil.copytree(source/'raw/early_synteny',root/'raw/prior_early_synteny'/data['review_id'])
    data['parent_review']=data['review_id'];data['review_id']=str(uuid4());data['created']=now()
    data['source_evidence_sha256']=sha256(source/'evidence.json')
    data.pop('early_synteny',None)
    for locus in data['loci']:locus.pop('early_synteny',None)
    error=attach(args,root,data)
    write_outputs(root,data)
    print('STOP — open '+str(root/'report.html'),file=sys.stderr)
    return 1 if error else 0
