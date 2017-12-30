from .web import CSSPrinter, _html_clsname

def _js_hover(s, script_tags=True):
    codestr = """
var classes = [%s];
var elms = {};
var n = {}, nclasses = classes.length;
function changeColor(classname, color) {
    var curN = n[classname];
    for(var i = 0; i < curN; i ++) {
        elms[classname][i].style.backgroundColor = color;
    }
}
for(var k = 0; k < nclasses; k ++) {
    var curClass = classes[k];
    elms[curClass] = document.getElementsByClassName(curClass);
    n[curClass] = elms[curClass].length;
    var curN = n[curClass];
    for(var i = 0; i < curN; i ++) {
        elms[curClass][i].onmouseover = function() {
            changeColor(this.className, "LightBlue");
        };
        elms[curClass][i].onmouseout = function() {
            changeColor(this.className, "inherit");
        };
    }
};
""" % s  # from https://stackoverflow.com/a/12786869/790973
    if script_tags:
        return '<script type="text/javascript">' + codestr + '</script>'
    else:
        return codestr

class JSPrinter(CSSPrinter):
    """ Prints javascript-enabled HTML representaions """

    def _print_Reaction(self, rxn, **kwargs):
        res = self._get_str_Reaction('html_name', 'html_arrow', substances, **kwargs)

    def _print(self, obj, **kwargs):
        return super(JSPrinter, self)._print(obj, **kwargs) + self._js()

    def _js(self, **kwargs):
        return _js_hover('"%s"' % '", "'.join(map(_html_clsname, self._get('substances', **kwargs))))


def javascript(obj, **settings):
    return JSPrinter(settings).doprint(obj)
