from itertools import chain


class Printer(object):
    _str = str  # set to unicdoe int UnicodePrinter for Python 2
    _default_settings = dict(
        with_param=True,
        fallback_print_fn=str,
        Reaction_param_separator='; ',
        Reaction_around_arrow=(' ', ' '),
        magnitude_fmt=lambda x: '%.3g' % x,
    )
    _default_setting_factories = dict(
        substances=dict,
        colors=dict,  # substance key -> (bg-color, border-color), 6 char hex colors
    )
    _default_setting_attrs = dict(
        Reaction_coeff_fmt='_str',
        Reaction_formula_fmt='_str',
        unit_fmt='_str',
    )
    printmethod_attr = None  # e.g. '_html' or '_unicode', allows object local printing logic

    def __init__(self, settings=None):
        self._settings = dict(self._default_settings, **(settings or {}))
        for k, v in self._default_setting_factories.items():
            if k not in self._settings:
                self._settings[k] = v()
        for k, v in self._default_setting_attrs.items():
            if k not in self._settings:
                self._settings[k] = getattr(self, v)
        for k in self._settings:
            if k not in chain(self._default_settings, self._default_setting_factories,
                              self._default_setting_attrs):
                raise ValueError("Unknown setting: %s (missing in default_settings)" % k)

    def _get(self, key, **kwargs):
        return kwargs.get(key, self._settings[key])

    def _print(self, obj, **kwargs):
        for cls in type(obj).__mro__:
            print_meth = '_print_' + cls.__name__
            if hasattr(self, print_meth):
                return getattr(self, print_meth)(obj, **kwargs)
            for PrintCls in self.__class__.__mro__:
                _attr = getattr(PrintCls, 'printmethod_attr', None)
                if _attr and hasattr(obj, _attr):
                    return getattr(obj, _attr)(self, **kwargs)
        fn = self._get('fallback_print_fn', **kwargs)
        if fn:
            return fn(obj)
        else:
            raise ValueError("Don't know how to print obj of type: %s" % type(obj))

    def doprint(self, obj):
        return self._print(obj)
