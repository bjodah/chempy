# -*- coding: utf-8 -*-
"""
Utilities for plotting with `bokeh <https://bokeh.pydata.org>`_.
"""
from __future__ import (absolute_import, division, print_function)

from collections import OrderedDict, defaultdict
from itertools import chain

from chempy.kinetics.ode import get_odesys
from chempy.units import to_unitless, linspace, logspace_from_lin


def integration_with_sliders(
        rsys, tend, c0, parameters, fig_kwargs=None, unit_registry=None, output_conc_unit=None,
        output_time_unit=None, slider_kwargs=None, x_axis_type="linear", y_axis_type="linear",
        integrate_kwargs=None, substitutions=None):

    import numpy as np
    from bokeh.plotting import Figure
    from bokeh.models import ColumnDataSource, HBox, VBoxForm
    from bokeh.models.widgets import Slider

    if slider_kwargs is None:
        slider_kwargs = {}
    odesys, state_keys, rarg_keys, p_units = get_odesys(
        rsys, unit_registry=unit_registry, substitutions=substitutions,
        output_conc_unit=output_conc_unit,
        output_time_unit=output_time_unit)[:4]
    if output_conc_unit is None:
        output_conc_unit = 1
    if output_time_unit is None:
        output_conc_unit = 1

    param_keys = list(chain(state_keys, rarg_keys))
    if x_axis_type == 'linear':
        tout = linspace(tend*0, tend)
    elif x_axis_type == 'log':
        tout = logspace_from_lin(tend*1e-9, tend)
    else:
        raise NotImplementedError("Unknown x_axis_type: %s" % x_axis_type)

    tout, Cout, info = odesys.integrate(tout, c0, parameters, **(integrate_kwargs or {}))
    sources = [ColumnDataSource(data={
        'tout': to_unitless(tout, output_time_unit),
        k: to_unitless(Cout[:, idx], output_conc_unit)
    }) for idx, k in enumerate(rsys.substances)]

    if fig_kwargs is None:
        Cmax = np.max(Cout)
        x_range = list(to_unitless([tend*0, tend], output_time_unit))
        y_range = list(to_unitless([Cmax*0, Cmax*1.1], output_conc_unit))
        fig_kwargs = dict(plot_height=400, plot_width=400, title="C vs t",
                          tools="crosshair,pan,reset,resize,save,wheel_zoom",
                          x_range=x_range, y_range=y_range, x_axis_type=x_axis_type,
                          y_axis_type=y_axis_type)
    plot = Figure(**fig_kwargs)

    colors = 'red green blue black cyan magenta'.split()
    for idx, k in enumerate(rsys.substances):
        plot.line('tout', k, source=sources[idx], line_width=3, line_alpha=0.6,
                  color=colors[idx % len(colors)])

    def _C(k):
        return to_unitless(c0[k], output_conc_unit)
    p_ul = [to_unitless(parameters[k], _u) for k, _u in zip(param_keys, p_units)]
    c0_widgets = OrderedDict([
        (k, Slider(
            title=k if output_conc_unit is 1 else k + ' / ' + output_conc_unit.dimensionality.unicode,
            value=_C(k), **slider_kwargs.get(k, dict(start=_C(k)/2, end=_C(k)*2, step=_C(k)/10))))
        for k in rsys.substances])

    def _dict_to_unitless(d, u):
        return {k: to_unitless(v, u) for k, v in d.items()}

    param_widgets = OrderedDict([
        (k, Slider(title=k if u is None else k + ' / ' + u.dimensionality.unicode,
                   value=v, **_dict_to_unitless(
                       slider_kwargs.get(k, dict(start=v/10, end=v*10, step=v/10)),
                       u)))
        for k, v, u in zip(param_keys, p_ul, p_units)])
    all_widgets = list(chain(param_widgets.values(), c0_widgets.values()))

    def update_data(attrname, old, new):
        _c0 = defaultdict(lambda: 0*output_conc_unit)
        for k, w in c0_widgets.items():
            _c0[k] = w.value * output_conc_unit
        _params = {}
        for (k, w), u in zip(param_widgets.items(), p_units):
            _params[k] = w.value if u is None else w.value * u
        _tout, _Cout, _info = odesys.integrate(tout, _c0, _params)
        for idx, k in enumerate(rsys.substances):
            sources[idx].data = {
                'tout': to_unitless(_tout, output_time_unit),
                k: to_unitless(_Cout[:, idx], output_conc_unit)
            }

    for w in all_widgets:
        w.on_change('value', update_data)

    inputs = VBoxForm(children=all_widgets)
    return HBox(children=[inputs, plot], width=800)
