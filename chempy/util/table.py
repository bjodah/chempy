# -*- coding: utf-8 -*-
"""
Convenience functions for presenting reaction systems in tables.
"""
from __future__ import (absolute_import, division, print_function)


import os
import shutil
import subprocess
import tempfile

from ..printing import latex
from ..kinetics.rates import RadiolyticBase
from ..units import to_unitless, get_derived_unit

tex_templates = {
    'document': {
        'default': r"""
\documentclass[a4paper,9pt]{article}
\pagestyle{empty}
\usepackage[paper=a4paper,margin=1cm]{geometry}
%(usepkg)s
\hypersetup{
  bookmarksnumbered=true,
  breaklinks=false,
  raiselinks=true,
  pdfborder={0 0 0},
  colorlinks=true,
  plainpages=false,
  pdfstartview={FitH},
  pdfcreator={LaTeX with hyperref package},
  citecolor=teal,
  linkcolor=red,
  urlcolor=blue,
}
\begin{document}
%(begins)s
%(table)s
%(ends)s
\end{document}
"""
    },
    'table': {
        'default': r"""
\begin{%(table_env)s}
\centering
\label{tab:%(label)s}
\caption[%(short_cap)s]{%(long_cap)s}
\begin{tabular}{%(alignment)s}
\toprule
%(header)s
\midrule
%(body)s
\bottomrule
\end{tabular}
\end{%(table_env)s}""",
        'longtable': r"""
\begin{%(table_env)s}{%(alignment)s}
\caption[%(short_cap)s]{%(long_cap)s
\label{tab:%(label)s}}\\
\toprule
%(header)s
\midrule
%(body)s
\bottomrule
\end{%(table_env)s}"""
    }
}


def render_tex_to_pdf(contents, texfname, pdffname, output_dir, save):
    """ Generates a pdf from a tex file by calling pdflatex

    Parameters
    ----------
    contents : str
    texfname : path
    pdffname : path
    output_dir : path
    save : path or bool or str(bool)

    """
    created_tempdir = False
    try:
        if output_dir is None:
            output_dir = tempfile.mkdtemp()
            created_tempdir = True
        texpath = os.path.join(output_dir, texfname)
        pdfpath = os.path.join(output_dir, pdffname)
        cmds = ['pdflatex', '-halt-on-error', '-interaction',
                'batchmode', texfname]
        with open(texpath, 'wt') as ofh:
            ofh.write(contents)
            ofh.flush()
        with open(pdfpath + '.out', 'wb') as logfile:
            p = subprocess.Popen(cmds, cwd=output_dir,
                                 stdout=logfile, stderr=logfile)
            retcode = p.wait()
            p = subprocess.Popen(cmds, cwd=output_dir,
                                 stdout=logfile, stderr=logfile)
            retcode += p.wait()
        if retcode:
            fmtstr = "{}\n returned with exit status {}"
            raise RuntimeError(fmtstr.format(' '.join(cmds), retcode))
        else:
            return pdfpath
    finally:
        if save is True or save == 'True':
            pass
        else:
            if save is False or save == 'False':
                if created_tempdir:
                    shutil.rmtree(output_dir)
            else:
                # interpret path to copy pdf to.
                if not os.path.samefile(pdfpath, save):
                    shutil.copy(pdfpath, save)


def rsys2tablines(rsys, rref0=1, coldelim=' & ',
                  tex=True, ref_fmt=None,
                  unit_registry=None, unit_fmt='{}', k_fmt='%.4g'):
    """
    Generates a table representation of a ReactionSystem.

    Parameters
    ----------
    rsys : ReactionSystem
    rref0 : integer
        default start of index counter (default: 1)
    coldelim : string
        column delimiter (default: ' & ')
    tex : bool
        use latex formated output (default: True)
    ref_fmt : string or callable
        format string of ``ref`` attribute of reactions
    unit_registry : unit registry
        optional (default: None)
    """
    if ref_fmt is None:
        def _doi(s):
            return r'\texttt{\href{http://dx.doi.org/'+s+'}{doi:'+s+'}}'

        def ref_fmt(s):
            if s is None:
                return 'None'
            if tex:
                if isinstance(s, dict):
                    return _doi(s['doi'])
                if s.startswith('doi:'):
                    return _doi(s[4:])
            return s

    def _wrap(s):
        if tex:
            return '\\ensuremath{' + s + '}'
        else:
            return s

    lines = []
    for ri, rxn in enumerate(rsys.rxns):
        rxn_ref = rxn.ref
        if isinstance(rxn.param, RadiolyticBase):
            if unit_registry is not None:
                kunit = get_derived_unit(unit_registry, 'radiolytic_yield')
                k = k_fmt % to_unitless(rxn.param.args[0], kunit)
                k_unit_str = (kunit.dimensionality.latex.strip('$') if tex
                              else kunit.dimensionality)
        else:
            if unit_registry is not None:
                kunit = (get_derived_unit(unit_registry,
                                          'concentration')**(1-rxn.order()) /
                         get_derived_unit(unit_registry, 'time'))
                try:
                    k = k_fmt % to_unitless(rxn.param, kunit)
                    k_unit_str = (kunit.dimensionality.latex.strip('$') if tex
                                  else kunit.dimensionality)
                except Exception:
                    k, k_unit_str = rxn.param.equation_as_string(k_fmt, tex)
            else:
                k_unit_str = '-'
                if isinstance(k_fmt, str):
                    k = k_fmt % rxn.param
                else:
                    k = k_fmt(rxn.param)
        latex_kw = dict(with_param=False)
        if tex:
            latex_kw['substances'] = rsys.substances
            latex_kw['Reaction_around_arrow'] = ('}}' + coldelim + '\\ensuremath{{',
                                                 '}}' + coldelim + '\\ensuremath{{')
        else:
            latex_kw['Reaction_around_arrow'] = (coldelim,)*2
            latex_kw['Reaction_arrow'] = '->'
        lines.append(coldelim.join([
            str(rref0+ri),
            ('\\ensuremath{%s}' if tex else '%s') % latex(rxn, **latex_kw),
            _wrap(k),
            unit_fmt.format(_wrap(k_unit_str)),
            ref_fmt(rxn_ref) if callable(ref_fmt) else ref_fmt.format(rxn_ref)
        ]))

    return lines


def rsys2table(rsys, table_template=None, table_template_dict=None,
               param_name='Rate constant', **kwargs):
    r"""
    Renders user provided table_template with table_template_dict which
    also has 'body' entry generated from `rsys2tablines`.

    Defaults is LaTeX table requiring booktabs package to be used
    (add \usepackage{booktabs} to preamble).

    Parameters
    ----------
    rsys : ReactionSystem
    table_template : string
    table_tempalte_dict : dict used to render table_template (excl. "body")
    param_name : str
        Column header for parameter column
    longtable : bool
        use longtable in defaults. (default: False)
    **kwargs :
        passed onto rsys2tablines

    """
    siunitx = kwargs.pop('siunitx', False)
    line_term = r' \\'
    defaults = {
        'table_env': 'longtable' if kwargs.pop(
            'longtable', False) else 'table',
        'alignment': 'llllSll' if siunitx else 'lllllll',
        'header': kwargs.get('coldelim', ' & ').join([
            'Id.', 'Reactants', '', 'Products', '{%s}' % param_name,
            'Unit', 'Ref'
        ]) + line_term,
        'short_cap': rsys.name,
        'long_cap': rsys.name,
        'label': (rsys.name or 'None').lower()
    }

    if table_template_dict is None:
        table_template_dict = defaults
    else:
        for k, v in defaults:
            if k not in table_template_dict:
                table_template_dict[k] = v

    if 'body' in table_template_dict:
        raise KeyError("There is already a 'body' key in table_template_dict")
    if 'k_fmt' not in kwargs:
        kwargs['k_fmt'] = r'\num{%.4g}' if siunitx else '%.4g'
    table_template_dict['body'] = (line_term + '\n').join(rsys2tablines(
        rsys, **kwargs)
    ) + line_term

    if table_template is None:
        if table_template_dict['table_env'] == 'longtable':
            table_template = tex_templates['table']['longtable']
        else:
            table_template = tex_templates['table']['default']

    return table_template % table_template_dict


def rsys2pdf_table(rsys, output_dir=None, doc_template=None,
                   doc_template_dict=None, save=True, landscape=False,
                   **kwargs):
    """
    Convenience function to render a ReactionSystem as
    e.g. a pdf using e.g. pdflatex.

    Parameters
    ----------
    rsys : ReactionSystem
    output_dir : path to output directory
        (default: system's temporary folder)
    doc_template : string
        LaTeX boiler plate temlpate including preamble,
        document environment etc.
    doc_template_dict : dict (string -> string)
        dict used to render temlpate (excl. 'table')
    longtable : bool
        use longtable in defaults. (default: False)
    **kwargs :
        passed on to `rsys2table`
    """
    if doc_template is None:
        doc_template = tex_templates['document']['default']
    lscape = ['pdflscape' if landscape == 'pdf' else 'lscape'] if landscape else []
    _pkgs = [
        'booktabs', 'amsmath', ('pdftex,colorlinks,unicode=True', 'hyperref')
    ] + lscape
    if kwargs.get('longtable', False):
        _pkgs += ['longtable']
    if kwargs.get('siunitx', False):
        _pkgs += ['siunitx']
    _envs = ['tiny'] + (['landscape'] if landscape else [])
    defaults = {
        'usepkg': '\n'.join([(r'\usepackage' + (
            '[%s]' if isinstance(pkg, tuple) else '') + '{%s}') % pkg for pkg in _pkgs]),
        'begins': '\n'.join([r'\begin{%s}' % env for env in _envs]),
        'ends': '\n'.join([r'\end{%s}' % env for env in _envs[::-1]])
    }

    if doc_template_dict is None:
        doc_template_dict = defaults
    else:
        for k, v in defaults:
            if k not in doc_template_dict:
                doc_template_dict[k] = v

    if 'table' in doc_template_dict:
        raise KeyError("There is already a 'table' key in doc_template_dict")
    doc_template_dict['table'] = rsys2table(rsys, **kwargs)

    contents = doc_template % doc_template_dict

    if isinstance(save, str) and save.endswith('.pdf'):
        texfname = save.rstrip('.pdf') + '.tex'
        pdffname = save
    else:
        texfname = 'output.tex'
        pdffname = 'output.pdf'
    return render_tex_to_pdf(contents, texfname, pdffname, output_dir, save)
