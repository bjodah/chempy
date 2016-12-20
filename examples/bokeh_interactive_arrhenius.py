"""
Interactive kinetics app with sliders.
Start by runing:
    $ bokeh serve interactive.py
Add --show argument or navigate to:
    http://localhost:5006/interactive
"""
from collections import defaultdict
import sys

from chempy import Reaction, ReactionSystem
from chempy.arrhenius import ArrheniusParam
from chempy.util.bkh import integration_with_sliders


def get_rsys(Af=1e16, Ab=1.5e15, Ea=72e3, Er=-12e3):
    fw = Reaction({'Fe+3', 'SCN-'}, {'FeSCN+2'}, ArrheniusParam(Af, Ea))
    bw = Reaction({'FeSCN+2'}, {'Fe+3', 'SCN-'}, ArrheniusParam(Ab, Ea-Er))
    return ReactionSystem([fw, bw], 'Fe+3 SCN- FeSCN+2'.split())


if __name__.startswith('bk_'):
    from bokeh.io import curdoc
    curdoc().add_root(integration_with_sliders(
        get_rsys(), tend=3,
        c0=defaultdict(float, {'Fe+3': 3e-3, 'SCN-': 1.5e-3, 'FeSCN+2': .1e-3}),
        parameters={'temperature': 298.15},
        slider_kwargs={'temperature': dict(start=273.15, end=313.15, step=.05)}
    ))
elif __name__ == '__main__':
    import warnings
    warnings.warn("Run using 'bokeh serve %s'" % __file__)
    sys.exit(1)
