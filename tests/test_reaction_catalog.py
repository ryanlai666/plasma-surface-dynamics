import importlib.util
import json
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    'reaction_check', ROOT / 'scripts/check_reaction_catalog.py'
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_catalog_conserves_atoms_and_charge():
    catalog = json.loads((ROOT / 'configs/reaction_candidates.json').read_text())
    assert module.validate(catalog) == 45
    assert all(not r['active_in_kmc'] and not r['rate_available'] for r in catalog['reactions'])


def test_rejects_missing_hydrogen_product():
    catalog = json.loads((ROOT / 'configs/reaction_candidates.json').read_text())
    del catalog['reactions'][0]['products']['H2']
    with pytest.raises(ValueError, match='Unbalanced'):
        module.validate(catalog)


def test_rejects_charge_imbalance_even_with_atom_balance():
    catalog = dict(
        species={
            'neutral': dict(composition={'Ar': 1}, charge=0),
            'ion': dict(composition={'Ar': 1}, charge=1),
        },
        reactions=[dict(id='bad_ionization', reactants={'neutral': 1}, products={'ion': 1})],
    )
    with pytest.raises(ValueError, match='Unbalanced'):
        module.validate(catalog)
