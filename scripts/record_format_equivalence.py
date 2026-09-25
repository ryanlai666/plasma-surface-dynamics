"""Record formatting-only code changes in the provenance ledger.

For every tracked Python file that differs from ``--base`` (a git revision),
this script checks that the old and new files have the same syntax
fingerprint.  Matching files are added to ``provenance/code_format_ledger.json``
so that result records carrying the old byte hash remain verifiable.  Files
whose logic changed are reported and removed from the ledger; any saved
results that record their old hash must be regenerated.

    python scripts/record_format_equivalence.py --base e0aba51
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from plasma_surface.provenance import LEDGER, sha256_bytes, syntax_fingerprint_source  # noqa: E402


def git(*args: str) -> bytes:
    return subprocess.check_output(["git", *args], cwd=ROOT)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--base", required=True, help="git revision holding the old bytes")
    args = parser.parse_args()
    changed = (
        git("diff", "--name-only", args.base, "--", "plasma_surface/*.py", "scripts/*.py")
        .decode()
        .split()
    )
    ledger = json.loads(LEDGER.read_text(encoding="utf-8")) if LEDGER.exists() else {}
    files = ledger.setdefault("files", {})
    rejected = []
    for name in changed:
        path = ROOT / name
        if not path.exists():
            continue
        try:
            old = git("show", f"{args.base}:{name}")
        except subprocess.CalledProcessError:
            continue  # new file: nothing recorded against it
        new = path.read_bytes()
        fingerprint = syntax_fingerprint_source(new)
        if syntax_fingerprint_source(old) != fingerprint:
            rejected.append(name)
            files.pop(name, None)  # no longer a formatting-only successor
            continue
        entry = files.setdefault(name, {"predecessor_sha256": []})
        old_sha = sha256_bytes(old)
        if old_sha not in entry["predecessor_sha256"]:
            entry["predecessor_sha256"].append(old_sha)
        entry["syntax_sha256"] = fingerprint
        entry["current_sha256"] = sha256_bytes(new)
    ledger["description"] = (
        "Formatting-only code changes. Each predecessor hash is a file version whose "
        "syntax fingerprint equals syntax_sha256; see plasma_surface/provenance.py."
    )
    ledger["files"] = dict(sorted(files.items()))
    LEDGER.parent.mkdir(exist_ok=True)
    LEDGER.write_text(json.dumps(ledger, indent=2) + "\n", encoding="utf-8")
    print(f"recorded {len(files)} files in {LEDGER.relative_to(ROOT)}")
    if rejected:
        print(
            "logic changed (removed from ledger; regenerate dependent results):",
            *rejected,
            sep="\n  ",
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
