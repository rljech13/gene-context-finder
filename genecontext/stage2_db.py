"""Explicit, offline-reusable protein catalogues for stage 2 identification."""
import gzip
import hashlib
import json
from pathlib import Path
import re
import sqlite3
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from .context_reference import reference_records
from .sequence import sha256
from .stage1 import dump, new_output, now
from .tools import Runner, tool_versions


def verify_db(path):
    root = Path(path).resolve()
    info = json.loads((root/'receipt.json').read_text())
    if info.get('schema') != 'stage2_protein_db_v1':
        raise ValueError('Not a prepared stage 2 protein database')
    if {p.name for p in root.glob('blastdb.*')} != {name for name in info['files'] if name.startswith('blastdb.')}:
        raise ValueError('Protein database file set changed')
    for name, digest in info['files'].items():
        p = (root/name).resolve()
        if root not in p.parents or not p.is_file() or sha256(p) != digest:
            raise ValueError('Protein database changed: '+name)
    return root, info


def prepare(args):
    root = new_output(args.out)
    runner = Runner(root)
    if not tool_versions(['makeblastdb'])['makeblastdb']['path']:
        raise ValueError('Preparing the protein database requires makeblastdb')
    download = None
    source = Path(args.fasta).resolve() if args.fasta else root/'download.fasta.gz'
    if not args.fasta:
        query = '(taxonomy_id:7742) AND (reviewed:true)'
        url = 'https://rest.uniprot.org/uniprotkb/stream?'+urlencode(dict(
            query=query, format='fasta', compressed='true', includeIsoform='true'))
        partial = source.with_suffix('.gz.part')
        with urlopen(Request(url, headers={'User-Agent':'gene-context-finder'}), timeout=120) as response:
            release = response.headers.get('X-UniProt-Release', '')
            if not release:
                raise ValueError('UniProt did not identify its release; incomplete download is preserved')
            download = dict(url=url, query=query, release=release,
                release_date=response.headers.get('X-UniProt-Release-Date', ''),
                total_results=response.headers.get('X-Total-Results', ''), include_isoforms=True)
            with partial.open('wb') as out:
                for block in iter(lambda:response.read(1024*1024), b''):
                    out.write(block)
        partial.rename(source)
    con = sqlite3.connect(root/'metadata.sqlite')
    con.execute('CREATE TABLE proteins(id TEXT PRIMARY KEY, original_id TEXT UNIQUE, gene TEXT, description TEXT, taxon TEXT, length INTEGER)')
    kept = skipped = 0
    with (root/'proteins.faa').open('w') as out:
        for original, header, seq in reference_records(source):
            if len(seq)<10 or not set(seq)<=set('ACDEFGHIKLMNPQRSTVWYXBZUO'):
                skipped += 1
                continue
            pid = 'R'+hashlib.sha256(original.encode()).hexdigest()[:24]
            name = re.search(r'\bGN=([^\s]+)', header)
            taxon = re.search(r'\bOX=(\d+)', header)
            try:
                con.execute('INSERT INTO proteins VALUES(?,?,?,?,?,?)',
                    (pid, original, name.group(1) if name else '', header, taxon.group(1) if taxon else '', len(seq)))
            except sqlite3.IntegrityError as error:
                raise ValueError('Duplicate reference identifier or hash collision: '+original) from error
            out.write('>'+pid+'\n'+seq+'\n'); kept += 1
    con.commit(); con.close()
    if not kept:
        raise ValueError('No usable proteins in reference FASTA')
    runner.run(['makeblastdb','-in',root/'proteins.faa','-dbtype','prot','-out',root/'blastdb'],root/'makeblastdb.stdout')
    files = [p for p in root.iterdir() if p.is_file() and (p.name.startswith('blastdb.') or p.name in {'proteins.faa','metadata.sqlite'})]
    info = dict(schema='stage2_protein_db_v1', created=now(), proteins=kept,
        skipped_invalid=skipped, source=str(source), source_sha256=sha256(source),
        source_scope='reviewed vertebrate UniProtKB + available isoforms' if download else 'user-supplied FASTA; taxonomy is not inferred',
        download=download, commands=runner.commands, files={p.name:sha256(p) for p in files})
    dump(root/'receipt.json',info)
    verify_db(root)
    print('Prepared protein database: '+str(root))
    return 0
