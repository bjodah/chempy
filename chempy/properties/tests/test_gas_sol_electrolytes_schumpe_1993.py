from chempy.units import default_units as u
from ..gas_sol_electrolytes_schumpe_1993 import lg_solubility_ratio


def test_lg_solubility_ratio():
    lgr = lg_solubility_ratio({'Br-': 0.05*u.molar, 'Na+': 0.050*u.molar}, 'N2O', units=u)
    assert lgr != 0   # TODO: calculate by hand the reference value
