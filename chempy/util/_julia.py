from collections import OrderedDict
from chempy import Reaction
from chempy.kinetics.rates import MassAction, RadiolyticBase
from chempy.units import to_unitless, default_units as u


def jl_dict(od):
    return "Dict([%s])" % ", ".join(['(:%s, %.4g)' % (k, v) for k, v in od.items()])


def _r(r, p, substmap, parmap, *, unit_conc, unit_time,
       variables=None):
    pk, = r.param.unique_keys
    if isinstance(r.param, MassAction):
        ratcoeff = to_unitless(
            p[pk], unit_conc**(1-r.order())/unit_time
        )
        if not r.inact_reac:
            r_str = '{}, {}'.format(parmap[pk], r.string(substances=substmap, with_param=False,
                                                         Reaction_arrow='-->', Reaction_coeff_space=''))
        else:
            all_keys = r.keys()
            reac_stoichs = r.all_reac_stoich(all_keys)
            act_stoichs = r.active_reac_stoich(all_keys)
            rate = '*'.join([parmap[pk]] + [('%s^%d' % (substmap[k], v)) if v > 1 else substmap[k]
                                            for k, v in zip(all_keys, act_stoichs) if v > 0])
            r2 = Reaction(dict([(k, v) for k, v in zip(all_keys, reac_stoichs) if v]), r.prod)
            r_str = '{}, {}'.format(rate, r2.string(substances=substmap, with_param=False,
                                                    Reaction_arrow='\u21D2', Reaction_coeff_space=''))
    elif isinstance(r.param, RadiolyticBase):
        ratcoeff = to_unitless(
            p[pk]*variables['doserate']*variables['density'],
            unit_conc/unit_time
        )
        assert not r.reac and not r.inact_reac and not r.inact_prod
        (prod, n), = r.prod.items()
        assert n == 1
        r_str = ('{}, 0 \u21D2 {}' if ratcoeff > 0 else '{}, {} \u21D2 0').format(
            parmap[pk], substmap[prod])
    else:
        raise NotImplementedError("Whats that?")
    return r_str, pk, abs(ratcoeff)


class DiffEqBioJl:
    _template_body = """\
{name} = @{crn_macro} begin
    {reactions}
end {parameters}
{post}
"""

    defaults = dict(unit_conc=u.molar, unit_time=u.second)

    def __init__(self, *, rxs, pars, substance_key_map, parmap, **kwargs):
        self.rxs = rxs
        self.pars = pars
        self.substance_key_map = substance_key_map
        self.parmap = parmap
        self.unit_conc = kwargs.get('unit_conc', self.defaults['unit_conc'])
        self.unit_time = kwargs.get('unit_time', self.defaults['unit_time'])
        self.latex_names = kwargs.get('latex_names', None)

    @classmethod
    def from_rsystem(cls, rsys, par_vals, *, variables=None, substance_key_map=lambda i, sk: 'y%d' % i, **kwargs):
        if not isinstance(substance_key_map, dict):
            substance_key_map = {sk: substance_key_map(si, sk) for si, sk in enumerate(rsys.substances)}
        parmap = dict([(r.param.unique_keys[0], 'p%d' % i) for i, r in enumerate(rsys.rxns)])
        rxs, pars = [], OrderedDict()
        for r in rsys.rxns:
            rs, pk, pv = _r(r, par_vals, substance_key_map, parmap, variables=variables,
                            unit_conc=kwargs.get('unit_conc', cls.defaults['unit_conc']),
                            unit_time=kwargs.get('unit_time', cls.defaults['unit_time']))
            rxs.append(rs)
            if pk in pars:
                raise ValueError("Are you sure (sometimes intentional)?")
            pars[parmap[pk]] = pv
        if 'latex_names' not in kwargs:
            kwargs['latex_names'] = {k: s.latex_name for k, s in rsys.substances.items()}
        return cls(rxs=rxs, pars=pars, substance_key_map=substance_key_map, parmap=parmap, **kwargs)

    def render_body(self, sparse_jac=False):
        name = 'rn'
        return self._template_body.format(
            crn_macro='min_reaction_network' if sparse_jac else 'reaction_network',
            name=name,
            reactions='\n    '.join(self.rxs),
            parameters=' '.join(self.pars),
            post="addodes!({}, sparse_jac=True)".format(name) if sparse_jac else ''
        )

    def render_setup(self, *, ics, atol, tex=True, tspan=None):
        export = ""
        export += "p = %s\n" % jl_dict(self.pars)
        export += "ics = %s\n" % jl_dict(OrderedDict({self.substance_key_map[k]: v for
                                                      k, v in to_unitless(ics, u.molar).items()}))
        if atol:
            export += "abstol_d = %s\n" % jl_dict({self.substance_key_map[k]: v for
                                                   k, v in to_unitless(atol, u.molar).items()})
            export += "abstol = Array([get(abstol_d, k, 1e-10) for k=keys(speciesmap(rn))])"
        if tex:
            export += "subst_tex = Dict([%s])\n" % ", ".join(
                '(:%s, ("%s", "%s"))' % (v, k, self.latex_names[k]) for
                k, v in self.substance_key_map.items())
        if tspan:
            export += """\
tspan = (0., %12.5g)
u0 = Array([get(ics, k, 1e-28) for k=keys(speciesmap(rn))])
parr = Array([p[k] for k=keys(paramsmap(rn))])
oprob = ODEProblem(rn, u0, tspan, parr)
"""
        return export

    def render_solve(self):
        return ("sol = solve(oprob, reltol=1e-9, abstol=abstol, Rodas5(), "
                "callback=PositiveDomain(ones(length(u0)), abstol=abstol))")
