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
from chempy.kinetics.rates import MassAction
from chempy.util.bkh import integration_with_sliders


def get_rsys():
    r1 = Reaction({'A'}, {'B'}, MassAction([4./100], unique_keys=['k1']), name='R1: A cons.')
    r2 = Reaction({'B', 'C'}, {'A', 'C'}, MassAction([1e4], unique_keys=['k2']), name='R2: A reform.')
    r3 = Reaction({'B': 2}, {'B', 'C'}, MassAction([3e7], unique_keys=['k3']), name='R3: C form.')
    return ReactionSystem([r1, r2, r3])


if __name__.startswith('bk_'):
    # TODO: better xlim, ylim, ranges for c0
    from bokeh.io import curdoc
    curdoc().add_root(integration_with_sliders(
        get_rsys(), tend=1e3,
        c0=defaultdict(float, {'A': 1}),
        parameters={'k1': 4./100, 'k2': 1e4, 'k3': 3e7},
        x_axis_type='log', y_axis_type='log',
        integrate_kwargs=dict(integrator='cvode')
    ))
elif __name__ == '__main__':
    import warnings
    warnings.warn("Run using 'bokeh serve %s'" % __file__)
    sys.exit(1)
