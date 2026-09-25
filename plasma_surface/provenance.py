"""Code provenance checks that survive formatting-only changes.

Result records store the SHA-256 of the code that produced them.  A byte hash
breaks when a file is merely reformatted, even though the program is
unchanged.  This module accepts a recorded hash when either

1. the current file bytes still match it, or
2. the formatting ledger lists the recorded hash as a predecessor of the
   current file, and the current file's syntax fingerprint equals the one
   stored in the ledger.

The syntax fingerprint is the SHA-256 of the parsed Python AST (without
line/column attributes and with docstring indentation normalised), so it is
identical for any two files that differ only in layout, quoting or comments.
Any change to logic, constants or names changes the fingerprint and the check
fails, which is the intended behaviour: a new program needs new results.

The ledger is written by ``scripts/record_format_equivalence.py``, which
verifies each old file (read from git) against its reformatted successor
before recording it.  The predecessors remain available in git history.
"""

from __future__ import annotations

import ast
import hashlib
import inspect
import json
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "provenance" / "code_format_ledger.json"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str | Path) -> str:
    return sha256_bytes(_resolve(path).read_bytes())


def syntax_fingerprint_source(source: str | bytes) -> str:
    """Hash of the Python AST, independent of layout, quotes and comments."""
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            body = node.body
            if (
                body
                and isinstance(body[0], ast.Expr)
                and isinstance(body[0].value, ast.Constant)
                and isinstance(body[0].value.value, str)
            ):
                body[0].value.value = inspect.cleandoc(body[0].value.value).strip()
    return sha256_bytes(ast.dump(tree, include_attributes=False).encode())


def syntax_fingerprint(path: str | Path) -> str:
    return syntax_fingerprint_source(_resolve(path).read_bytes())


@lru_cache(maxsize=1)
def load_ledger() -> dict:
    if not LEDGER.exists():
        return {}
    return json.loads(LEDGER.read_text(encoding="utf-8")).get("files", {})


def matches_recorded(path: str | Path, recorded_sha256: str) -> bool:
    """True if ``path`` has the recorded hash, or is a verified reformat of that code."""
    resolved = _resolve(path)
    if sha256_bytes(resolved.read_bytes()) == recorded_sha256:
        return True
    entry = load_ledger().get(_relative(resolved))
    if not entry or recorded_sha256 not in entry["predecessor_sha256"]:
        return False
    return syntax_fingerprint(resolved) == entry["syntax_sha256"]


def _resolve(path: str | Path) -> Path:
    path = Path(str(path).replace("\\", "/"))
    return path if path.is_absolute() else ROOT / path


def _relative(path: Path) -> str | None:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return None  # outside the repository, so never in the ledger
