import json

from plasma_surface.provenance import (
    LEDGER,
    ROOT,
    matches_recorded,
    sha256_bytes,
    syntax_fingerprint,
    syntax_fingerprint_source,
)


def test_fingerprint_ignores_layout_but_not_logic():
    compact = "def f(x):\n    '''Doc.'''\n    a=x+1;return a*2 # comment\n"
    formatted = 'def f(x):\n    """Doc."""\n    a = x + 1\n    return a * 2\n'
    changed = "def f(x):\n    '''Doc.'''\n    a=x+1;return a*3\n"
    assert syntax_fingerprint_source(compact) == syntax_fingerprint_source(formatted)
    assert syntax_fingerprint_source(compact) != syntax_fingerprint_source(changed)


def test_byte_hash_still_accepted_for_unledgered_files(tmp_path):
    p = tmp_path / 'artifact.txt'
    p.write_bytes(b'data')
    assert matches_recorded(p, sha256_bytes(b'data'))
    assert not matches_recorded(p, sha256_bytes(b'other'))


def test_every_ledger_entry_matches_the_current_code():
    files = json.loads(LEDGER.read_text(encoding='utf-8'))['files']
    assert files
    for name, entry in files.items():
        assert syntax_fingerprint(ROOT / name) == entry['syntax_sha256'], name
        for old in entry['predecessor_sha256']:
            assert matches_recorded(name, old)
