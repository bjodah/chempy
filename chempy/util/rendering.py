import re
from .parsing import get_parsing_context

class TemplateEvaluator:
    def __init__(self, pattern=r'\${(.*?)}', fmt='${%s}', globals_=None):
        self.mark = re.compile(pattern)
        self.fmt = fmt
        if globals_ is None:
            globals_ = get_parsing_context()
        self.globals_ = globals_
        
    def __call__(self, template, **kwargs):
        for item in self.mark.findall(template):
            ev = eval(item, dict(self.globals_, **kwargs))
            template=template.replace(self.fmt % item, str(ev))
        return template

eval_template = TemplateEvaluator()
