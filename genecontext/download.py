"""Download pinned genome inputs with official MD5 verification and provenance."""
from __future__ import annotations

import csv
import hashlib
from http.client import HTTPException
import json
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import URLError
from urllib.parse import unquote, urlsplit
from urllib.request import Request, urlopen

from .sequence import write_tsv


def md5(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checksum_entries(text):
    entries = {}
    for line in text.splitlines():
        parts = line.split(maxsplit=1)
        if len(parts) != 2 or not re.fullmatch(r"[0-9a-fA-F]{32}", parts[0]):
            continue
        # NCBI checksum lists include nested assembly_structure files with
        # repeated basenames (e.g. chr2acc). Keep their relative paths distinct.
        filename = parts[1].removeprefix("*").removeprefix("./")
        if filename in entries and entries[filename] != parts[0].lower():
            raise ValueError(f"Conflicting official MD5 entries: {filename}")
        entries[filename] = parts[0].lower()
    if not entries:
        raise ValueError("Official checksum file contains no MD5 entries")
    return entries


def transfer(url, destination, expected, retries=3):
    """Resume a partial file when supported, verify bytes, then rename atomically."""
    destination = Path(destination)
    if destination.exists():
        if md5(destination) == expected:
            return "reused_verified"
        raise ValueError(f"Existing file fails official MD5; preserved unchanged: {destination}. "
                         "Use a new output directory or inspect/remove that file explicitly.")
    partial = destination.with_name(destination.name + ".part")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if partial.exists() and md5(partial) == expected:
        partial.replace(destination)
        return "completed_verified_partial"
    for attempt in range(retries):
        offset = partial.stat().st_size if partial.exists() else 0
        headers = {"User-Agent": "sequence-gene-context/0.1.0"}
        if offset:
            headers["Range"] = f"bytes={offset}-"
        try:
            with urlopen(Request(url, headers=headers), timeout=60) as response:
                resume = offset > 0 and response.status == 206
                if resume and not response.headers.get("Content-Range", "").startswith(f"bytes {offset}-"):
                    raise ValueError(f"Incorrect HTTP range response: {url}")
                # Servers ignoring Range return 200: restart instead of duplicating bytes.
                with open(partial, "ab" if resume else "wb") as handle:
                    for block in iter(lambda: response.read(1024 * 1024), b""):
                        handle.write(block)
            observed = md5(partial)
            if observed != expected:
                raise ValueError(f"Official MD5 mismatch for {url}: expected {expected}, got {observed}; "
                                 f"unverified bytes remain at {partial}")
            partial.replace(destination)
            return "downloaded_verified"
        except (OSError, URLError, HTTPException) as error:
            if attempt + 1 == retries:
                raise RuntimeError(f"Download failed after {retries} attempts: {url}: {error}") from error
            time.sleep(min(2**attempt, 4))
    raise RuntimeError(f"Download did not finish: {url}")


def load_sources(path, only=None):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"genome_id", "genome_url"}.issubset(reader.fieldnames or []):
            raise ValueError("Download TSV requires genome_id and genome_url columns")
        rows = [{key: value or "" for key, value in row.items()} for row in reader if any(row.values())]
    seen = set()
    for row in rows:
        gid = row["genome_id"]
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", gid) or gid in seen:
            raise ValueError(f"Unsafe or duplicate genome_id: {gid!r}")
        seen.add(gid)
        if not row["genome_url"]:
            raise ValueError(f"{gid}: missing genome_url")
        for key in ("genome_url", "gff_url", "protein_url", "cds_url", "md5_url"):
            url = row.get(key, "")
            if url and (urlsplit(url).scheme not in {"http", "https"} or not urlsplit(url).netloc):
                raise ValueError(f"{gid}: invalid HTTP(S) URL in {key}")
    if only:
        unknown = set(only) - {r["genome_id"] for r in rows} - {r.get("accession", "") for r in rows}
        if unknown:
            raise ValueError("Unknown --only identifiers/accessions: " + ", ".join(sorted(unknown)))
        rows = [r for r in rows if r["genome_id"] in only or r.get("accession") in only]
    if not rows:
        raise ValueError("No genomes selected for download")
    return rows


def run_download(args):
    sources = load_sources(args.manifest, args.only)
    kinds = ["genome", "gff"] + (["protein", "cds"] if args.include_reference else [])
    plan = [dict(genome_id=row["genome_id"], accession=row.get("accession", ""),
                 files={kind: row[kind + "_url"] for kind in kinds if row.get(kind + "_url")},
                 annotation_available=bool(row.get("gff_url"))) for row in sources]
    if args.dry_run:
        print(json.dumps(plan, indent=2))
        return 0
    root = Path(args.out).resolve()
    root.mkdir(parents=True, exist_ok=True)
    ready, receipts, failures, checksums = [], [], [], {}
    suffixes = dict(genome="genome.fna", gff="annotation.gff3", protein="reference.protein.faa", cds="reference.cds.fna")
    for row in sources:
        gid, downloaded = row["genome_id"], {}
        receipt = dict(genome_id=gid, accession=row.get("accession", ""), files=[])
        receipts.append(receipt)
        directory = root / gid
        directory.mkdir(exist_ok=True)
        try:
            for kind in kinds:
                url = row.get(kind + "_url", "")
                if not url:
                    continue
                checksum_url = row.get("md5_url") or url.rsplit("/", 1)[0] + "/md5checksums.txt"
                if checksum_url not in checksums:
                    with urlopen(Request(checksum_url, headers={"User-Agent": "sequence-gene-context/0.1.0"}), timeout=60) as response:
                        raw = response.read(16 * 1024 * 1024 + 1)
                    if len(raw) > 16 * 1024 * 1024:
                        raise ValueError("Checksum file unexpectedly exceeds 16 MiB")
                    text = raw.decode("utf-8")
                    checksums[checksum_url] = (text, checksum_entries(text))
                text, entries = checksums[checksum_url]
                filename = Path(unquote(urlsplit(url).path)).name
                if filename not in entries:
                    raise ValueError(f"No official MD5 for {filename} in {checksum_url}")
                checksum_name = "source_" + hashlib.sha256(checksum_url.encode()).hexdigest()[:12] + ".md5"
                (directory / checksum_name).write_text(text)
                target = directory / (suffixes[kind] + (".gz" if filename.endswith(".gz") else ""))
                print(f"[{gid}] {kind}: downloading/verifying {url}", file=sys.stderr, flush=True)
                status = transfer(url, target, entries[filename])
                receipt["files"].append(dict(kind=kind, url=url, path=str(target),
                    checksum_url=checksum_url, md5=entries[filename], bytes=target.stat().st_size, status=status))
                downloaded[kind] = str(target.relative_to(root))
            ready.append(dict(genome_id=gid, fasta=downloaded["genome"], gff3=downloaded.get("gff", "")))
            receipt["status"] = "ready"
        except (OSError, ValueError, RuntimeError, URLError, HTTPException) as error:
            receipt.update(status="error", error=str(error))
            failures.append(dict(genome_id=gid, error=str(error)))
            print(f"[{gid}] ERROR: {error}; continuing", file=sys.stderr)
    write_tsv(root / "genomes.tsv", ready, ["genome_id", "fasta", "gff3"])
    (root / "download_receipt.json").write_text(json.dumps(dict(
        checked_utc=datetime.now(timezone.utc).isoformat(), source_manifest=str(Path(args.manifest).resolve()),
        genomes=receipts, errors=failures, ready=len(ready)), indent=2) + "\n")
    print(f"{len(ready)} ready genomes; {len(failures)} failures. Input manifest: {root/'genomes.tsv'}", file=sys.stderr)
    return 1 if failures else 0
