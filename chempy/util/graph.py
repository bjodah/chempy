# -*- coding: utf-8 -*-
"""
Convenince functions for representing reaction systems as graphs.
"""
import os
import subprocess
import shutil
import tempfile


def rsys2dot(rsys, tex=False, rprefix='r', rref0=1,
             nodeparams='[label="{}" shape=diamond]', colors=('maroon', 'darkgreen')):
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
    lines = ['digraph ' + str(rsys.name) + '{\n']
    ind = '  '  # indentation

    def add_vertex(key, num, reac):
        snum = str(num) if num > 1 else ''
        name = getattr(rsys.substances[key], 'latex_name' if tex else 'name')
        lines.append(ind + '"{}" -> "{}" [label ="{}",color={},fontcolor={}];\n'.format(
            *((name, rid, snum, colors[0], colors[0]) if reac else
              (rid, name, snum, colors[1], colors[1]))
        ))

    all_reac_stoichs = rsys.all_reac_stoichs()
    all_prod_stoichs = rsys.all_prod_stoichs()

    for ri, rxn in enumerate(rsys.rxns):
        rid = rprefix + str(ri+rref0)
        lines.append(ind + '{')
        lines.append(ind*2 + 'node ' + nodeparams.format(rxn.name or rid))
        lines.append(ind*2 + rid)
        lines.append(ind + '}\n')
        for idx, key in enumerate(rsys.substances):
            num = all_reac_stoichs[ri, idx]
            if num == 0:
                continue
            add_vertex(key, num, True)
        for idx, key in enumerate(rsys.substances):
            num = all_prod_stoichs[ri, idx]
            if num == 0:
                continue
            add_vertex(key, num, False)
    lines.append('}\n')
    return lines


def rsys2graph(rsys, fname, output_dir=None, prog=None, save=False, **kwargs):
    """
    Convenience function to call `rsys2dot` and write output to file
    and render the graph

    Parameters
    ----------
    rsys: ReactionSystem
    fname: str
        filename
    output_dir: str (optional)
        path to directory (default: temporary directory)
    prog: str (optional)
        default: 'dot'
    save: bool
        removes temporary directory if False, default: False
    \*\*kwargs:
        parameters to pass along to `rsys2dot`

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
