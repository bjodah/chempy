"""
Interactive kinetics app with sliders.
Start by runing:
    $ bokeh serve interactive.py
Add --show argument or navigate to:
    http://localhost:5006/interactive
"""
from collections import defaultdict

from bokeh.io import curdoc

from chempy import Reaction, ReactionSystem
from chempy.kinetics.rates import MassAction
from chempy.util.bkh import integration_with_sliders


def get_rsys(kf, kb):
    rf = MassAction([kf], unique_keys=['kf'])
    rb = MassAction([kb], unique_keys=['kb'])
    fw = Reaction({'Fe+3', 'SCN-'}, {'FeSCN+2'}, rf)
    bw = Reaction({'FeSCN+2'}, {'Fe+3', 'SCN-'}, rb)
    return ReactionSystem([fw, bw], 'Fe+3 SCN- FeSCN+2'.split())


if __name__.startswith('bk_'):
    curdoc().add_root(integration_with_sliders(
        get_rsys(3, .3), tend=3,
        c0=defaultdict(float, {'Fe+3': .9, 'SCN-': .7}),
        parameters={'kf': 3, 'kb': .3}
    ))
else:
    import warnings
    warnings.warn("Run using 'bokeh serve %s'" % __file__)
