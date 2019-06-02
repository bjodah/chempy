from collections import OrderedDict
from chempy import Reaction
from chempy.kinetics.rates import MassAction, RadiolyticBase
from chempy.units import to_unitless, default_units as u


_diffeqbiojl_tmplt = """\
{name} = @reaction_network begin
    {reactions}
end {parameters}
"""


def _r(r, p, doserate, substmap, parmap, density=998*u.kg/u.m3, unit_conc=u.molar, unit_time=u.second):
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
            p[pk]*doserate*density,
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


def to_diffeqbiojl(arm, arm_extra, *, doserate):
    sbstmap = dict(zip(arm.substances, (chr(i) for i in range(ord('A'), ord('A')+999))))
    parmap = dict(zip([r.param.unique_keys[0] for r in arm.rxns], ('k%d' % i for i in range(1, 999))))
    rxs, pars = [], OrderedDict()
    for r in arm.rxns:
        rs, pk, pv = _r(r, arm_extra['params'], doserate, sbstmap, parmap)
        rxs.append(rs)
        if pk in pars:
            raise ValueError("Are you sure (sometimes intentional)?")
        pars[parmap[pk]] = pv

    return dict(
        rn=_diffeqbiojl_tmplt.format(
            name='rn',
            reactions='\n    '.join(rxs),
            parameters=' '.join(pars)
        ),
        sbstmap=sbstmap,
        parmap=parmap,
        pars=pars
    )


def export2julia(armor_rsys, armor_extra, *, ics, kw2=None):
    debj = to_diffeqbiojl(armor_rsys, armor_extra, **(kw2 or {}))
    export = debj['rn']

    def p_odj(od):
        return "Dict([%s])" % ", ".join(['(:%s, %.4g)' % (k, v) for k, v in od.items()])

    # str(od).replace('Ordered', '').replace("',", ',').replace("'", ":")
    export += "p = %s\n" % p_odj(debj['pars'])
    export += "ics = %s\n" % p_odj(OrderedDict({debj['sbstmap'][k]: v for
                                                k, v in to_unitless(ics, u.molar).items()}))
    export += "subst_tex = Dict([%s])\n" % ", ".join(
        '(:%s, ("%s", "%s"))' % (v, k, armor_rsys.substances[k].latex_name) for
        k, v in debj['sbstmap'].items())
    return export
