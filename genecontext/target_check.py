"""Small, offline protein-reference checks to support manual target review."""
from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import re

from .context_reference import reference_records
from .review import new_root, page
from .sequence import sha256, union_length, write_fasta, write_tsv
from .comparison import table
from .tools import Runner, tool_versions


FIELDS = ['candidate_id', 'rank', 'reference_id', 'reference_name', 'family',
          'candidate_aa', 'reference_aa', 'identity', 'candidate_coverage',
          'reference_coverage', 'bitscore', 'evalue', 'hsps_used']


def proteins(path, prefix, *, reference=False):
    sequences, metadata, seen = [], {}, set()
    for identifier, header, sequence in reference_records(path):
        if identifier in seen:
            raise ValueError('Duplicate protein identifier: '+identifier)
        seen.add(identifier)
        stops = sequence.count('*')
        if not sequence or not set(sequence) <= set('ACDEFGHIKLMNPQRSTVWYXBZJUO*'):
            raise ValueError('Invalid protein sequence: '+identifier)
        if reference and stops:
            raise ValueError('Reference protein contains internal stops: '+identifier)
        key = prefix+f'{len(sequences)+1:07d}'
        sequences.append((key, sequence.replace('*', 'X')))
        name = re.search(r'\bGN=(\S+)', header)
        metadata[key] = dict(id=identifier, description=header, protein_aa=len(sequence),
                             internal_stops=stops, name=name.group(1) if name else identifier)
    return sequences, metadata


def labels(path, reference_metadata):
    if not path:
        return {}
    expected = {r['id'] for r in reference_metadata.values()}
    result = {}
    for row in table(path):
        identifier, family = (row.get('reference') or '').strip(), (row.get('family') or '').strip()
        if identifier not in expected or identifier in result or not family:
            raise ValueError('Reference labels require unique known reference IDs and nonempty families: '+identifier)
        result[identifier] = family
    if set(result) != expected:
        raise ValueError('Reference labels must cover every reference protein')
    return result


def summarize(hsps, candidate_length, reference_length):
    """Greedy, nonoverlapping local HSPs; not a global protein alignment."""
    query, subject, chosen = [], [], []
    for h in sorted(hsps, key=lambda r: -r['bitscore']):
        q, s = (h['qstart'], h['qend']), (h['sstart'], h['send'])
        if any(max(q[0], a) < min(q[1], b) for a, b in query):
            continue
        if any(max(s[0], a) < min(s[1], b) for a, b in subject):
            continue
        query.append(q); subject.append(s); chosen.append(h)
    length = sum(h['length'] for h in chosen)
    return dict(identity=sum(h['identity']*h['length'] for h in chosen)/length,
                candidate_coverage=union_length(query)/candidate_length,
                reference_coverage=union_length(subject)/reference_length,
                bitscore=sum(h['bitscore'] for h in chosen),
                evalue=min(h['evalue'] for h in chosen), hsps_used=len(chosen))


def check_targets(args):
    from . import __version__
    candidates, cm = proteins(args.proteins, 'Q')
    references, rm = proteins(args.references, 'R', reference=True)
    families = labels(args.reference_labels, rm)
    root = Path(args.out).resolve()
    inputs = [Path(args.proteins).resolve(), Path(args.references).resolve()]
    if args.reference_labels:
        inputs.append(Path(args.reference_labels).resolve())
    if any(root == p or root in p.parents for p in inputs):
        raise ValueError('Output must not contain the input files')
    versions = tool_versions(['blastp'])
    if not versions['blastp']['path']:
        raise ValueError('check-targets requires BLAST+ blastp on PATH')
    root = new_root(root, inputs[0])
    write_fasta(root/'candidates.faa', candidates)
    write_fasta(root/'references.faa', references)
    for filename, metadata in [('candidate_ids.tsv', cm), ('reference_ids.tsv', rm)]:
        rows = [dict(search_id=k, **v, family=families.get(v['id'], '')) for k, v in metadata.items()]
        write_tsv(root/filename, rows, ['search_id', 'id', 'name', 'protein_aa', 'internal_stops', 'family', 'description'])
    runner = Runner(root)
    try:
        runner.run(['blastp', '-query', root/'candidates.faa', '-subject', root/'references.faa',
                    '-num_threads', args.threads, '-evalue', '1e-5', '-max_target_seqs', len(references),
                    '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen'],
                   root/'blastp.tsv')
        grouped = defaultdict(list)
        for line in (root/'blastp.tsv').read_text().splitlines():
            f = line.split('\t')
            if len(f) != 12 or f[0] not in cm or f[1] not in rm:
                raise ValueError('Unexpected protein-reference BLASTP output')
            qlen, slen = cm[f[0]]['protein_aa'], rm[f[1]]['protein_aa']
            qa, qb, sa, sb = map(int, [f[4], f[5], f[6], f[7]])
            if int(f[10]) != qlen or int(f[11]) != slen or not (1 <= qa <= qb <= qlen and 1 <= sa <= sb <= slen):
                raise ValueError('BLASTP coordinates/lengths do not match input proteins')
            grouped[f[0], f[1]].append(dict(identity=float(f[2])/100, length=int(f[3]),
                qstart=qa-1, qend=qb, sstart=sa-1, send=sb, evalue=float(f[8]), bitscore=float(f[9])))
        by_candidate = defaultdict(list)
        for (qid, rid), hsps in grouped.items():
            ref = rm[rid]
            by_candidate[qid].append(dict(candidate_id=cm[qid]['id'], reference_id=ref['id'],
                reference_name=ref['name'], family=families.get(ref['id'], ''),
                candidate_aa=cm[qid]['protein_aa'], reference_aa=ref['protein_aa'],
                **summarize(hsps, cm[qid]['protein_aa'], ref['protein_aa'])))
        hits, summary = [], []
        for qid, candidate in cm.items():
            matches = sorted(by_candidate[qid], key=lambda r: (-r['bitscore'], r['evalue'], r['reference_id']))
            for rank, hit in enumerate(matches, 1):
                hit['rank'] = rank
                hits.append(hit)
            best = matches[0] if matches else {}
            alternative = next((h for h in matches[1:] if h['family'] != best['family']), None) if families else None
            ratio = alternative['bitscore']/best['bitscore'] if alternative else ''
            flags = []
            if candidate['internal_stops']:
                flags.append('internal_stop')
            if not best:
                flags.append('no_reference_hit')
            elif min(best['candidate_coverage'], best['reference_coverage']) < .5:
                flags.append('limited_coverage')
            if ratio != '' and ratio >= .95:
                flags.append('close_other_family')
            summary.append(dict(candidate_id=candidate['id'], candidate_aa=candidate['protein_aa'],
                flags=';'.join(flags) or 'reference_match',
                **{'best_'+key: best.get(key, '') for key in ['reference_id', 'reference_name', 'family', 'identity',
                                                             'candidate_coverage', 'reference_coverage', 'bitscore']},
                alternative_family=alternative['family'] if alternative else '',
                alternative_reference=alternative['reference_id'] if alternative else '',
                alternative_score_ratio=ratio))
        write_tsv(root/'hits.tsv', hits, FIELDS)
        write_tsv(root/'summary.tsv', summary, list(summary[0]))
        page(root/'report.html', 'Protein reference check — evidence for manual review',
             '<p>Local BLASTP matches against the supplied reference panel. Coverage and identity are fractions (0–1). '
             'An internal stop is searched as X and remains flagged. A no-hit result does not establish absence. '
             'Family labels are supplied by the user; score ratios are not orthology probabilities. '
             'References outside this panel have not been tested. No review decisions are changed.</p>'
             '<p><a href="hits.tsv">All reference matches</a> · <a href="summary.tsv">Summary TSV</a> · '
             '<a href="reference_ids.tsv">Reference panel and labels</a> · <a href="manifest.json">Provenance</a></p>'
             '<p>Flags: limited_coverage = either coverage below 0.5; close_other_family = another labelled family '
             'scores at least 0.95 of the best. These are inspection aids, not acceptance criteria.</p>',
             summary, list(summary[0]))
    except (OSError, ValueError, RuntimeError) as error:
        (root/'manifest.json').write_text(json.dumps(dict(version=__version__, status='failed', error=str(error),
            commands=runner.commands), indent=2)+'\n')
        raise
    manifest = dict(version=__version__, status='completed', finished=datetime.now(timezone.utc).isoformat(),
        inputs={str(p): sha256(p) for p in inputs}, tools=versions, commands=runner.commands,
        candidates=len(cm), references=len(rm), matching_candidates=len(by_candidate),
        parameters=dict(evalue=1e-5, limited_coverage_below=.5, close_other_family_ratio=.95,
                        aggregation='greedy HSPs nonoverlapping on both axes'),
        limits='Reference-panel similarity evidence only; no automatic gene naming, orthology or review approval.')
    manifest['outputs_sha256'] = {str(p.relative_to(root)): sha256(p) for p in root.iterdir() if p.is_file()}
    (root/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return 0
