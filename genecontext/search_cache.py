"""Shared, verified search inputs; independent of query panels and reports."""
import json
from pathlib import Path
import shutil
import tempfile
import time

from .early_synteny import cache_lock, seal_cache, token, verified_cache
from .sequence import sha256
from .tools import Genome, tool_versions


def cache_root(args):
    """Serial and per-genome workers must resolve to the same shared root."""
    return Path(getattr(args, 'neighbor_cache', None) or
                Path(args.out).resolve().parent/'_genecontext_cache'/'early_synteny').resolve()


def fingerprint(path):
    s = Path(path).stat()
    return s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns


def tool_identity(name, versions=None):
    value = (versions or tool_versions([name]))[name]
    if not value.get('path'):
        raise ValueError('Missing search-cache tool: '+name)
    return dict(path=str(Path(value['path']).resolve()), version=value['version'])


def prepared(root, kind, identity, runner, build):
    """Publish a complete cache atomically; never replace corrupted evidence."""
    started = time.perf_counter()
    key = token(identity)
    destination = Path(root)/'search_inputs'/kind/key
    destination.parent.mkdir(parents=True, exist_ok=True)
    with cache_lock(destination.with_suffix('.lock')):
        reused = verified_cache(destination, identity)
        if reused:
            expected = set(json.loads((destination/'cache.json').read_text())['outputs'])
            actual = {str(p.relative_to(destination)) for p in destination.rglob('*')
                      if p.is_file() and p.name != 'cache.json'}
            if actual != expected:
                raise ValueError('Unexpected shared search-cache files: '+str(destination))
        if not reused:
            if destination.exists():
                raise ValueError('Incomplete shared search cache: '+str(destination))
            temporary = Path(tempfile.mkdtemp(prefix=key+'.building-', dir=destination.parent))
            commands_start = len(runner.commands)
            try:
                build(temporary)
                seal_cache(temporary, identity)
                temporary.rename(destination)
                # Logs move with the cache. Receipts must point to existing files
                # and retain the command that was actually run in the temp dir.
                for command in runner.commands[commands_start:]:
                    for field in ['output', 'stderr']:
                        value = command.get(field)
                        if value and Path(value).is_relative_to(temporary):
                            command[field] = str(destination/Path(value).relative_to(temporary))
            except BaseException:
                if temporary.exists():
                    shutil.rmtree(temporary)
                raise
    return destination, dict(key=key, path=str(destination), reused=reused,
                             elapsed_seconds=time.perf_counter()-started,
                             identity=identity, verification='sha256_all_outputs')


def prepare_genome(source, root, runner, source_sha256=None, versions=None):
    """Own one expanded FASTA/faidx; supplied digest must be caller-verified."""
    source = Path(source).resolve()
    before = fingerprint(source)
    digest = source_sha256 or sha256(source)
    identity = dict(schema='gcf_genome_index_v1', source_sha256=digest,
                    samtools=tool_identity('samtools', versions),
                    implementation_sha256=sha256(Path(__file__)),
                    genome_reader_sha256=sha256(Path(__file__).with_name('tools.py')))

    def build(folder):
        genome = Genome(source, folder, runner)
        try:
            # Ordinary FASTA previously used a symlink. The shared input owns
            # its bytes so later edits to the source cannot alter old indices.
            if genome.path.is_symlink():
                genome.path.unlink()
                shutil.copyfile(source, genome.path)
        finally:
            genome.close()
        if fingerprint(source) != before or sha256(source) != digest:
            raise ValueError('Source FASTA changed during index preparation')

    folder, receipt = prepared(root, 'genomes', identity, runner, build)
    if fingerprint(source) != before:
        raise ValueError('Source FASTA changed while checking shared index')
    receipt['fasta_sha256'] = json.loads((folder/'cache.json').read_text())['outputs']['genome.fa']
    return Genome.from_index(folder/'genome.fa'), receipt


def prepare_blast_db(fasta, root, runner, dbtype, source_sha256=None, versions=None):
    """Reuse exact inputs; a supplied digest must be from a verified read."""
    if dbtype not in {'prot', 'nucl'}:
        raise ValueError('Unsupported BLAST database type: '+str(dbtype))
    fasta = Path(fasta).resolve()
    before = fingerprint(fasta)
    digest = source_sha256 or sha256(fasta)
    identity = dict(schema='gcf_blast_database_v1', source_sha256=digest, dbtype=dbtype,
                    makeblastdb=tool_identity('makeblastdb', versions),
                    implementation_sha256=sha256(Path(__file__)))

    def build(folder):
        runner.run(['makeblastdb', '-in', fasta, '-dbtype', dbtype,
                    '-out', folder/'db'], folder/'makeblastdb.stdout')
        required = ['pin', 'phr', 'psq'] if dbtype == 'prot' else ['nin', 'nhr', 'nsq']
        # BLAST may split large databases into numbered volumes.
        for extension in required:
            if not list(folder.glob('db*.'+extension)):
                raise ValueError('Incomplete BLAST database: '+extension)
        if fingerprint(fasta) != before or sha256(fasta) != digest:
            raise ValueError('FASTA changed during BLAST database preparation')

    folder, receipt = prepared(root, 'blast_'+dbtype, identity, runner, build)
    if fingerprint(fasta) != before:
        raise ValueError('FASTA changed while checking BLAST database')
    return folder/'db', receipt
