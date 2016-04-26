"""
Interactive kinetics app with sliders (with units).
Start by runing:
    $ bokeh serve interactive.py
Add --show argument or navigate to:
    http://localhost:5006/interactive
"""
from collections import defaultdict

from bokeh.io import curdoc
from chempy.util.bkh import integration_with_sliders
from chempy.units import SI_base_registry, default_units as u

from bokeh_interactive import get_rsys


if __name__.startswith('bk_'):
    kf, kb = 3/u.molar/u.s, .3/u.s
    curdoc().add_root(integration_with_sliders(
        get_rsys(kf, kb), tend=3*u.s,
        c0=defaultdict(lambda: 0*u.molar, {'Fe+3': .9*u.molar, 'SCN-': .7*u.molar}),
        parameters={'kf': kf, 'kb': kb},
        unit_registry=SI_base_registry,
        output_conc_unit=u.molar,
        output_time_unit=u.second
    ))
else:
    import warnings
    warnings.warn("Run using 'bokeh serve %s'" % __file__)
