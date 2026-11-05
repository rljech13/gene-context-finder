"""A reusable, checksum-pinned Swiss-Prot fallback for unannotated windows."""
from datetime import datetime, timezone
from http.client import HTTPException
import json
from pathlib import Path
import re
import sys
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

from .download import md5, transfer
from .sequence import open_text, sha256


SWISSPROT_BASE = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/"
DEFAULT_CACHE = str(Path.home() / ".cache" / "gene-context" / "references")


def swissprot(cache, base_url=SWISSPROT_BASE, allow_download=True):
    """Reuse the pinned local release offline; fetch only when the cache is empty.

    base_url is injectable for local HTTP integration tests, not a CLI option.
    A new cache directory explicitly opts into a newer upstream release.
    """
    root = Path(cache).expanduser().resolve() / "swissprot"
    root.mkdir(parents=True, exist_ok=True)
    receipt = root / "receipt.json"
    if receipt.exists():
        info = json.loads(receipt.read_text())
        if not isinstance(info, dict) or not {"bytes", "md5", "sha256", "release"}.issubset(info):
            raise ValueError(f"Invalid Swiss-Prot cache receipt; preserved unchanged: {receipt}. Use a new --reference-cache.")
        path = root / "uniprot_sprot.fasta.gz"
        if (not path.is_file() or path.stat().st_size != info["bytes"] or
                md5(path) != info["md5"] or sha256(path) != info["sha256"]):
            raise ValueError(f"Swiss-Prot cache fails verification; preserved unchanged: {root}. Use a new --reference-cache.")
        return path, dict(info, cache_status="reused_verified", path=str(path))
    if not allow_download:
        raise ValueError("Swiss-Prot is not cached. Run prepare-context or setup-predictor --with-reference before --offline, or use --context-reference none.")
    # Prevent simultaneous first downloads from writing the same .part file.
    lock = root / "download.lock"
    try:
        handle = lock.open("x")
    except FileExistsError as error:
        raise RuntimeError(f"Swiss-Prot cache is being prepared: {lock}. Retry after that run finishes; remove a stale lock only after checking it.") from error
    try:
        with handle:
            handle.write(datetime.now(timezone.utc).isoformat())
        url = base_url + "RELEASE.metalink"
        with urlopen(Request(url, headers={"User-Agent": "gene-context-finder"}), timeout=60) as response:
            raw = response.read(1024 * 1024 + 1)
        if len(raw) > 1024 * 1024:
            raise ValueError("Swiss-Prot release metadata exceeds 1 MiB")
        document = ET.fromstring(raw)
        ns = {"m": "http://www.metalinker.org/"}
        entries = [item for item in document.findall("m:files/m:file", ns)
                   if item.get("name") == "uniprot_sprot.fasta.gz"]
        if len(entries) != 1:
            raise ValueError("Swiss-Prot release metadata has no unique FASTA entry")
        entry = entries[0]
        expected = entry.findtext("m:verification/m:hash[@type='md5']", namespaces=ns) or ""
        size = int(entry.findtext("m:size", namespaces=ns) or "0")
        release = document.findtext("m:version", namespaces=ns) or ""
        if not re.fullmatch(r"[a-fA-F0-9]{32}", expected) or size <= 0 or not release:
            raise ValueError("Incomplete Swiss-Prot release metadata")
        destination = root / "uniprot_sprot.fasta.gz"
        print(f"[context] Downloading Swiss-Prot {release}: {size/1e6:.1f} MB; cache: {root}", file=sys.stderr, flush=True)
        status = transfer(base_url + destination.name, destination, expected.lower())
        if destination.stat().st_size != size:
            raise ValueError("Swiss-Prot file size differs from official release metadata")
        (root / "RELEASE.metalink").write_bytes(raw)
        info = dict(database="UniProtKB/Swiss-Prot", release=release, url=base_url + destination.name,
                    metadata_url=url, md5=expected.lower(), sha256=sha256(destination), bytes=size,
                    downloaded_utc=datetime.now(timezone.utc).isoformat(), license="CC BY 4.0",
                    attribution="UniProt Consortium", license_url="https://www.uniprot.org/help/license")
        temporary = receipt.with_suffix(".json.part")
        temporary.write_text(json.dumps(info, indent=2) + "\n")
        temporary.replace(receipt)
        return destination, dict(info, cache_status=status, path=str(destination))
    except (HTTPException, ET.ParseError) as error:
        raise RuntimeError(f"Could not prepare Swiss-Prot: {error}") from error
    finally:
        lock.unlink()


def reference_records(path):
    """Yield identifier, description and sequence without materialising the file."""
    header, sequence = None, []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header.split()[0], header, "".join(sequence).upper().removesuffix("*")
                header, sequence = line[1:].strip(), []
                if not header:
                    raise ValueError(f"Empty context protein header: {path}")
            elif line.strip():
                if header is None:
                    raise ValueError(f"Sequence before context protein header: {path}")
                sequence.append("".join(line.split()))
    if header is None:
        raise ValueError(f"Empty context protein FASTA: {path}")
    yield header.split()[0], header, "".join(sequence).upper().removesuffix("*")


def add_references(path, proteins, metadata, evidence, skip_invalid=False):
    kept, skipped, seen = 0, 0, set()
    for identifier, header, sequence in reference_records(path):
        if identifier in seen:
            raise ValueError(f"Duplicate context protein identifier: {identifier}")
        seen.add(identifier)
        if len(sequence) < 30 or not set(sequence) <= set("ACDEFGHIKLMNPQRSTVWYX"):
            if skip_invalid:
                skipped += 1
                continue
            raise ValueError(f"Invalid context reference protein: {identifier}")
        pid = f"R{len(metadata)+1:07d}"
        match = re.search(r"\bGN=([^\s]+)", header)
        proteins[pid] = sequence
        metadata[pid] = dict(genome_id="external_reference", gene_id=identifier,
                             name=match.group(1) if match else identifier, evidence=evidence,
                             description=header)
        kept += 1
    if not kept:
        raise ValueError(f"No usable context reference proteins: {path}")
    return dict(usable_proteins=kept, skipped_short_or_nonstandard=skipped)
