"""
Interactive kinetics app with sliders (with units).
Start by runing:
    $ bokeh serve interactive.py
Add --show argument or navigate to:
    http://localhost:5006/interactive
"""
from collections import defaultdict
import sys

from chempy.util.bkh import integration_with_sliders
from chempy.units import SI_base_registry, default_units as u

from bokeh_interactive_arrhenius import get_rsys


if __name__.startswith('bk_'):
    from bokeh.io import curdoc
    Af, Ab, Ea, Er = 1e16/u.molar/u.s, 1.5e15/u.s, 72e3*u.J/u.mol, -12e3*u.J/u.mol
    curdoc().add_root(integration_with_sliders(
        get_rsys(Af, Ab, Ea, Er), tend=3*u.s,
        c0=defaultdict(lambda: 0*u.molar, {'Fe+3': 3e-3*u.molar, 'SCN-': 1.5e-3*u.molar}),
        parameters={'temperature': 298.15*u.K},
        slider_kwargs={'temperature': dict(start=273.15*u.K, end=313.15*u.K, step=.05*u.K)},
        get_odesys_kw=dict(
            unit_registry=SI_base_registry,
            output_conc_unit=u.molar,
            output_time_unit=u.second
        )
    ))
elif __name__ == '__main__':
    import warnings
    warnings.warn("Run using 'bokeh serve %s'" % __file__)
    sys.exit(1)
