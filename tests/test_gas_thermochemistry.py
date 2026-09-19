import math
import pytest
from plasma_surface.gas_thermochemistry import gas_thermo, R_J_MOL_K, KJ_MOL_PER_EV


def test_source_table_values():
    h=gas_thermo('H2',400,100000)
    assert h['Cp_J_mol_K']==pytest.approx(29.18,abs=.01)
    assert h['S_standard_J_mol_K']==pytest.approx(139.2,abs=.1)
    w=gas_thermo('H2O',500,100000)
    assert w['Cp_J_mol_K']==pytest.approx(35.22,abs=.01)
    assert w['H_minus_H298_kJ_mol']==pytest.approx(6.92,abs=.01)


@pytest.mark.parametrize('species',['HF','SiF4','NH3','H2','HCl','H2O'])
def test_thermodynamic_derivatives_and_pressure(species):
    T=600.;d=.001
    center=gas_thermo(species,T,100000)
    plus=gas_thermo(species,T+d,100000)
    minus=gas_thermo(species,T-d,100000)
    assert (plus['H_minus_H298_kJ_mol']-minus['H_minus_H298_kJ_mol'])/(2*d)==pytest.approx(center['Cp_J_mol_K']/1000,rel=1e-7)
    low=gas_thermo(species,T,100000/math.e)
    assert center['mu_minus_H298_eV']-low['mu_minus_H298_eV']==pytest.approx(R_J_MOL_K*T/1000/KJ_MOL_PER_EV)
    assert not center['qualified_for_surface_rate']


def test_source_fit_limits_are_enforced():
    with pytest.raises(ValueError):gas_thermo('H2O',450,100000)
    for p in (0,-1,float('nan')):
        with pytest.raises(ValueError):gas_thermo('HF',450,p)
