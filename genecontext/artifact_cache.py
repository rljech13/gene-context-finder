"""Atomic, content-verified artifacts for expensive resumable steps."""
import json
from pathlib import Path
import shutil
import tempfile
import time

from .early_synteny import cache_lock, seal_cache, token, verified_cache


def cached(root, identity, build):
    started = time.perf_counter()
    key = token(identity)
    destination = Path(root) / key
    destination.parent.mkdir(parents=True, exist_ok=True)
    with cache_lock(destination.with_suffix('.lock')):
        reused = verified_cache(destination, identity)
        if reused:
            expected = set(json.loads((destination/'cache.json').read_text())['outputs'])
            actual = {str(p.relative_to(destination)) for p in destination.rglob('*')
                      if p.is_file() and p.name != 'cache.json'}
            if actual != expected:
                raise ValueError('Unexpected artifact-cache files: '+str(destination))
        else:
            if destination.exists():
                raise ValueError('Incomplete artifact cache; preserved unchanged: '+str(destination))
            temporary = Path(tempfile.mkdtemp(prefix=key+'.building-', dir=destination.parent))
            try:
                build(temporary)
                seal_cache(temporary, identity)
                temporary.rename(destination)
            finally:
                if temporary.exists():
                    shutil.rmtree(temporary)
    return destination, dict(key=key, path=str(destination), reused=reused,
                             elapsed_seconds=time.perf_counter()-started,
                             verification='sha256_all_outputs')
