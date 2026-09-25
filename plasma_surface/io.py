import csv
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path


def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("Cannot write empty table")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def manifest(config):
    source = Path(__file__).parent
    result = dict(
        config=config,
        python=platform.python_version(),
        packages={
            m: importlib.metadata.version(m)
            for m in ["numpy", "scipy", "scikit-learn", "matplotlib"]
        },
        source_sha256={
            f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(source.glob("*.py"))
        },
        data_origin="synthetic_uncalibrated_model",
    )
    if config.get("backend") == "cpp":
        from .native import library_path

        binary = library_path()
        build = binary.parent / "build_manifest.json"
        result["native_binary_sha256"] = hashlib.sha256(binary.read_bytes()).hexdigest()
        result["native_build"] = (
            json.loads(build.read_text(encoding="utf-8")) if build.exists() else None
        )
    return result
