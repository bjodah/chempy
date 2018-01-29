# -*- coding: utf-8 -*-
"""
Convenince functions for representing reaction systems as graphs.
"""
import os
import subprocess
import shutil
import tempfile


def rsys2dot(rsys, tex=False, rprefix='r', rref0=1, nodeparams='[label="{}",shape=diamond]',
             colors=('maroon', 'darkgreen'), penwidths=None, include_inactive=True):
    """
    Returns list of lines of DOT (graph description language)
    formated graph.

    Parameters
    ==========
    rsys: ReactionSystem
    tex: bool (default False)
        If set True, output will be LaTeX formated
    (Substance need to have latex_name attribute set)
    rprefix: string
        Reaction enumeration prefix, default: r
    rref0: integer
        Reaction enumeration inital counter value, default: 1
    nodeparams: string
        DOT formated param list, default: [label={} shape=diamond]

    Returns
    =======
    list of lines of DOT representation of the graph representation.

    """
    lines = ['digraph "' + str(rsys.name) + '" {\n']
    ind = '  '  # indentation
    if penwidths is None:
        penwidths = [1.0]*rsys.nr

    categories = rsys.categorize_substances(checks=())

    def add_substance(key):
        fc = 'black'
        if key in categories['depleted']:
            fc = colors[0]
        if key in categories['accumulated']:
            fc = colors[1]
        label = ('$%s$' if tex else '%s') % getattr(rsys.substances[key], 'latex_name' if tex else 'name')
        lines.append(ind + '"{key}" [fontcolor={fc} label="{lbl}"];\n'.format(key=key, fc=fc, lbl=label))

    for sk in rsys.substances:
        add_substance(sk)

    def add_vertex(key, num, reac, penwidth):
        snum = str(num) if num > 1 else ''
        fmt = ','.join(
            ['label="{}"'.format(snum)] +
            (['penwidth={}'.format(penwidth)] if penwidth != 1 else [])
        )
        lines.append(ind + '"{}" -> "{}" [color={},fontcolor={},{}];\n'.format(
            *((key, rid, colors[0], colors[0], fmt) if reac else
              (rid, key, colors[1], colors[1], fmt))
        ))

    if include_inactive:
        reac_stoichs = rsys.all_reac_stoichs()
        prod_stoichs = rsys.all_prod_stoichs()
    else:
        reac_stoichs = rsys.active_reac_stoichs()
        prod_stoichs = rsys.active_prod_stoichs()

    for ri, rxn in enumerate(rsys.rxns):
        rid = rprefix + str(ri+rref0)
        lines.append(ind + '{')
        lines.append(ind*2 + 'node ' + nodeparams.format(rxn.name or rid))
        lines.append(ind*2 + rid)
        lines.append(ind + '}\n')
        for idx, key in enumerate(rsys.substances):
            num = reac_stoichs[ri, idx]
            if num == 0:
                continue
            add_vertex(key, num, True, penwidths[ri])
        for idx, key in enumerate(rsys.substances):
            num = prod_stoichs[ri, idx]
            if num == 0:
                continue
            add_vertex(key, num, False, penwidths[ri])
    lines.append('}\n')
    return lines


def rsys2graph(rsys, fname, output_dir=None, prog=None, save=False, **kwargs):
    """
    Convenience function to call `rsys2dot` and write output to file
    and render the graph

    Parameters
    ----------
    rsys : ReactionSystem
    fname : str
        filename
    output_dir : str (optional)
        path to directory (default: temporary directory)
    prog : str (optional)
        default: 'dot'
    save : bool
        removes temporary directory if False, default: False
    \\*\\*kwargs :
        Keyword arguments passed along to py:func:`rsys2dot`.

    Returns
    -------
    str
        Outpath

    Examples
    --------
    >>> rsys2graph(rsys, sbstncs, '/tmp/out.png')  # doctest: +SKIP

    """

    lines = rsys2dot(rsys, **kwargs)
    created_tempdir = False
    try:
        if output_dir is None:
            output_dir = tempfile.mkdtemp()
            created_tempdir = True
        basename, ext = os.path.splitext(os.path.basename(fname))
        outpath = os.path.join(output_dir, fname)
        dotpath = os.path.join(output_dir, basename + '.dot')
        with open(dotpath, 'wt') as ofh:
            ofh.writelines(lines)
        if ext == '.tex':
            cmds = [prog or 'dot2tex']
        else:
            cmds = [prog or 'dot', '-T'+outpath.split('.')[-1]]
        p = subprocess.Popen(cmds + [dotpath, '-o', outpath])
        retcode = p.wait()
        if retcode:
            fmtstr = "{}\n returned with exit status {}"
            raise RuntimeError(fmtstr.format(' '.join(cmds), retcode))
        return outpath
    finally:
        if save is True or save == 'True':
            pass
        else:
            if save is False or save == 'False':
                if created_tempdir:
                    shutil.rmtree(output_dir)
            else:
                # interpret save as path to copy pdf to.
                shutil.copy(outpath, save)
