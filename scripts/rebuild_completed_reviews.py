"""Host-side presentation rebuild. No computation service or source checkout changes."""
from pathlib import Path
import json,hashlib,sys,zipfile
B=Path('/mnt/data/i.rubinshteyn/gene_context_experiments')
R=B/'presentation_methods/legacy_modern_review_20261120'
CODE=R/'code-v2'
if not CODE.exists():
    CODE.mkdir();zipfile.ZipFile(R/'code-v2.zip').extractall(CODE)
sys.path.insert(0,str(CODE))
from genecontext.legacy_review import build,package,sha
active=B/'results/pluripotency_queries_top20_20261111/delivery/interactive_review'
published=json.loads((active/'published.json').read_text())
plan=json.loads((B/'stage4_batches/outputs_20261119/plan.json').read_text())
completed={g['gene']:g for g in plan['genes']}
DEST=B/'presentations/legacy_modern_review_20261120'
DEST.mkdir(exist_ok=True);(DEST/'reports').mkdir(exist_ok=True)
entries=[]
for entry in published['reports']:
    if entry['stage']!=3:continue
    rid=entry['id'];source=(active/'reports'/rid).resolve();out=DEST/'reports'/rid
    if out.exists():
        print('exists',rid,flush=True);continue
    prior=completed.get(entry['gene']);confirmed=None
    if prior and prior['review_id']==entry['review_id']:
        confirmed=Path(prior['confirmed_decisions'])
        assert sha(confirmed.read_bytes())==prior['confirmed_decisions_sha256']
    record=build(source,out,entry['gene'],confirmed)
    record.update(id=rid,gene=entry['gene'],source_directory=str(source),directory=str(out),
                  old_public_target=str(Path(entry['public']).resolve()),old_archive_target=str(Path(entry['archive_alias']).resolve()),
                  public=entry['public'],archive_alias=entry['archive_alias'],url=entry['url'])
    (out/'build_record.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'id':rid,'models':record['original_models'],'confirmed':bool(confirmed)}),flush=True)
print('Builds complete',flush=True)
