"""Acquire a small versioned public DFT archive and record its provenance."""

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.request
import zipfile
from .io import write_csv, write_json

RECORD = "10211009"
EXPECTED_MD5 = "34644b3c55b44441dc64492759cb6d31"


def fetch_hcl(root="data/raw/si_hcl"):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    api = f"https://zenodo.org/api/records/{RECORD}"
    with urllib.request.urlopen(api, timeout=60) as response:
        metadata = json.load(response)
    entry = next(f for f in metadata["files"] if f["key"] == "data.zip")
    if entry["checksum"] != f"md5:{EXPECTED_MD5}" or entry["size"] > 1_000_000:
        raise ValueError("Public archive changed; review provenance before importing")
    if metadata["metadata"]["license"]["id"] != "cc-by-4.0":
        raise ValueError("Dataset license changed; review before importing")
    with urllib.request.urlopen(entry["links"]["self"], timeout=60) as response:
        payload = response.read(1_000_001)
    if hashlib.md5(payload).hexdigest() != EXPECTED_MD5:
        raise ValueError("Archive checksum mismatch")
    archive = root / "data.zip"
    archive.write_bytes(payload)
    write_json(root / "record.json", metadata)
    write_json(
        root / "provenance.json",
        dict(
            source=api,
            doi=metadata["metadata"]["doi"],
            license="CC-BY-4.0",
            creators=metadata["metadata"]["creators"],
            retrieved_utc=datetime.now(timezone.utc).isoformat(),
            sha256=hashlib.sha256(payload).hexdigest(),
            md5=EXPECTED_MD5,
            use="DFT geometry inventory only; no reaction barriers or etch rates inferred",
        ),
    )
    return inventory(archive, root / "structure_inventory.csv")


def inventory(archive, output):
    """Read archive members without extracting paths or executing third-party files."""
    rows = []
    with zipfile.ZipFile(archive) as z:
        if sum(f.file_size for f in z.infolist()) > 20_000_000:
            raise ValueError("Archive exceeds inventory size budget")
        for member in z.infolist():
            if (
                member.is_dir()
                or "__MACOSX" in member.filename
                or member.filename.split("/")[-1].startswith(".")
            ):
                continue
            content = z.read(member).decode("utf-8", errors="replace")
            lines = content.splitlines()
            species = counts = ""
            # VASP5 species and count lines; non-CONTCAR files remain in inventory.
            if len(lines) >= 8:
                try:
                    numbers = [int(n) for n in lines[6].split()]
                    names = lines[5].split()
                    if len(numbers) == len(names) and all(n >= 0 for n in numbers):
                        species, counts = " ".join(names), " ".join(map(str, numbers))
                except ValueError:
                    pass
            rows.append(
                dict(
                    member=member.filename,
                    bytes=member.file_size,
                    species=species,
                    counts=counts,
                    sha256=hashlib.sha256(z.read(member)).hexdigest(),
                    source_doi="10.5281/zenodo.10211009",
                )
            )
    write_csv(output, rows)
    return dict(
        files=len(rows), parsed_structures=sum(bool(r["species"]) for r in rows), output=str(output)
    )
