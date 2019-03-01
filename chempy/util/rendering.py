import re
from .parsing import get_parsing_context
from ..units import fold_constants


class TemplateEvaluator:
    def __init__(self, pattern=r'\${(.*?)}', fmt='${%s}', globals_=None, post_procs=()):
        self.mark = re.compile(pattern)
        self.fmt = fmt
        if globals_ is None:
            globals_ = get_parsing_context()
        self.globals_ = globals_
        self._post_procs = post_procs

    def _post_proc(self, arg):
        for pp in self._post_procs:
            arg = pp(arg)
        return arg

    def __call__(self, template, **kwargs):
        for item in self.mark.findall(template):
            ev = self._post_proc(eval(item, dict(self.globals_, **kwargs)))
            template = template.replace(self.fmt % item, str(ev))
        return template

eval_template = TemplateEvaluator(post_procs=(fold_constants, lambda x: str(x).replace(' ', '*')))
