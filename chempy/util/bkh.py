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
        rsys, tend, c0, parameters, fig_kwargs=None, slider_kwargs=None, conc_bounds=None,
        x_axis_type="linear", y_axis_type="linear", integrate_kwargs=None, odesys_extra=None,
        get_odesys_kw=None, integrate=None):
    """
    Parameters
    ----------
    rsys : ReactionSystem
    tend : float like
    c0 : dict
        Initial concentrations.
    parameters : dict
        Parameter values.
    fig_kwargs : dict
        Keyword-arguments passed to bokeh's ``Figure``.
    slider_kwargs : dict
        Keyword-arguments passed to bokeh's ``Slider``.
    conc_bounds : dict of dicts
        Mapping substance key to dict of bounds ('start', 'end', 'step').
    x_axis_type : str
    y_axis_type : str
    integrate_kwargs : dict
        Keyword-arguments passed to integrate.
    odesys_extra : tuple
        If odesys & extra have already been generated (avoids call to ``get_odesys``).
    get_odesys_kw : dict
        Keyword-arguments passed to ``get_odesys``.
    integrate : callback
        Defaults to ``odesys.integrate``.

    """

    import numpy as np
    from bokeh.plotting import Figure
    from bokeh.models import ColumnDataSource, Column, Row
    from bokeh.models.widgets import Slider

    if slider_kwargs is None:
        slider_kwargs = {}
    if get_odesys_kw is None:
        get_odesys_kw = {}
    if odesys_extra is None:
        odesys, extra = get_odesys(rsys, **get_odesys_kw)
    else:
        odesys, extra = odesys_extra
    if integrate is None:
        integrate = odesys.integrate

    state_keys, rarg_keys, p_units = [extra[k] for k in ('param_keys', 'unique', 'p_units')]
    output_conc_unit = get_odesys_kw.get('output_conc_unit', None)
    output_time_unit = get_odesys_kw.get('output_time_unit', None)
    unit_registry = get_odesys_kw.get('unit_registry', None)
    if output_conc_unit is None:
        if unit_registry is not None:
            raise ValueError("if unit_registry is given, output_conc_unit must also be given")
        output_conc_unit = 1
    if output_time_unit is None:
        if unit_registry is not None:
            raise ValueError("if unit_registry is given, output_time_unit must also be given")
        output_conc_unit = 1

    param_keys = list(chain(state_keys, rarg_keys))
    if x_axis_type == 'linear':
        tout = linspace(tend*0, tend)
    elif x_axis_type == 'log':
        tout = logspace_from_lin(tend*1e-9, tend)
    else:
        raise NotImplementedError("Unknown x_axis_type: %s" % x_axis_type)

    result = integrate(tout, c0, parameters, **(integrate_kwargs or {}))
    sources = [ColumnDataSource(data={
        'tout': to_unitless(result.xout, output_time_unit),
        k: to_unitless(result.yout[:, idx], output_conc_unit)
    }) for idx, k in enumerate(rsys.substances)]
    if fig_kwargs is None:
        Cmax = np.max(result.yout)
        x_range = list(to_unitless([result.xout[0], result.xout[-1]], output_time_unit))
        y_range = list(to_unitless([Cmax*0, Cmax*1.1], output_conc_unit))
        fig_kwargs = dict(plot_height=400, plot_width=400, title="C vs t",
                          tools="crosshair,pan,reset,save,wheel_zoom",
                          x_range=x_range, y_range=y_range, x_axis_type=x_axis_type,
                          y_axis_type=y_axis_type)
    plot = Figure(**fig_kwargs)

    colors = 'red green blue black cyan magenta'.split()
    for idx, k in enumerate(rsys.substances):
        plot.line('tout', k, source=sources[idx], line_width=3, line_alpha=0.6,
                  color=colors[idx % len(colors)])

    def _C(k):
        return to_unitless(c0[k], output_conc_unit)
    if p_units is None:
        p_units = [None]*len(param_keys)
    p_ul = [to_unitless(parameters[k], _u) for k, _u in zip(param_keys, p_units)]

    def _dict_to_unitless(d, u):
        return {k: to_unitless(v, u) for k, v in d.items()}

    c0_widgets = OrderedDict()
    for k in rsys.substances:
        if conc_bounds is not None and k in conc_bounds:
            if k in slider_kwargs:
                raise ValueError("Key '%s' both in slider_kwargs and conc_bounds" % k)
            slider_defaults = _dict_to_unitless(conc_bounds[k], output_conc_unit)
        else:
            ck = _C(k)
            if ck == 0:
                max_ = max(*[_C(k) for k in rsys.substances])
                slider_defaults = dict(start=0, end=max_, step=max_/100)
            else:
                slider_defaults = dict(start=_C(k)/2, end=_C(k)*2, step=_C(k)/10)
        c0_widgets[k] = Slider(
            title=k if output_conc_unit is 1 else k + ' / ' + output_conc_unit.dimensionality.unicode,
            value=_C(k), **slider_kwargs.get(k, slider_defaults)
        )

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
        _result = integrate(tout, _c0, _params)
        for idx, k in enumerate(rsys.substances):
            sources[idx].data = {
                'tout': to_unitless(_result.xout, output_time_unit),
                k: to_unitless(_result.yout[:, idx], output_conc_unit)
            }

    for w in all_widgets:
        w.on_change('value', update_data)

    inputs = Column(children=all_widgets)
    return Row(children=[inputs, plot], width=800)
