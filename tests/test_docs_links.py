"""Every relative link in the project's Markdown documents must resolve."""

import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
LINK = re.compile(r'\]\(([^)\s]+)\)')


def markdown_files():
    try:
        out = subprocess.check_output(['git', 'ls-files', '*.md'], cwd=ROOT, text=True)
        files = out.split()
    except (OSError, subprocess.CalledProcessError):
        files = [str(p.relative_to(ROOT)) for p in ROOT.rglob('*.md')]
    return sorted(f for f in files if not f.startswith(('third_party/', 'archive/', '.')))


@pytest.mark.parametrize('name', markdown_files())
def test_relative_links_resolve(name):
    path = ROOT / name
    broken = []
    for target in LINK.findall(path.read_text(encoding='utf-8', errors='replace')):
        if target.startswith(('http://', 'https://', 'mailto:', '#')):
            continue
        if not (path.parent / target.split('#')[0]).exists():
            broken.append(target)
    assert not broken, f'{name}: {broken}'
