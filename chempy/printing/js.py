import json
from .web import CSSPrinter, _html_clsname

_js_rsys_template = """
var cls_names_substances = %(cls_names_substances)s;
var substance_row_cls_irrel = %(substance_row_cls_irrel)s;
var elms = {};
var n = {}, nsubstances = cls_names_substances.length;
var nirrel = {};
function changeColor(classname, color) {
    var curN = n[classname];
    for(var i = 0; i < curN; i++) {
        elms[classname][i].style.backgroundColor = color;
    }
}
function toggleVisibility(classname_substance) {
    var curN = nirrel[classname_substance];
    for (var i=0; i<curN; i++) {
        var objs = document.getElementsByClassName(substance_row_cls_irrel[classname_substance][i]);
        for (var j=0; j<objs.length; ++j){
            objs[j].style.display = objs[j].style.display == "none" ? "table-row" : "none";
        }
    }
}
function resetTab(tab){
    tab.style.border = "0px";
    var rows = tab.getElementsByTagName('tr');
    [].forEach.call(rows, function(row){
        row.style.display = "table-row";
    });
    tab.getElementsByTagName('th')[0].innerHTML = tab.ori_header +
        "<br>(click on species to show a subset of reactions)";
};

for(var k = 0; k < nsubstances; k++) {
    var curClass = cls_names_substances[k];
    var curIrrel = substance_row_cls_irrel[k];
    elms[curClass] = document.getElementsByClassName(curClass);
    n[curClass] = elms[curClass].length;
    nirrel[curClass] = substance_row_cls_irrel[curClass].length;
    var curN = n[curClass];
    for(var i = 0; i < curN; i++) {
        elms[curClass][i].onmouseover = function() {
            changeColor(this.className, "LightBlue");
        };
        elms[curClass][i].onmouseout = function() {
            changeColor(this.className, "inherit");
        };
        elms[curClass][i].onclick = function() {
            var tab = this.closest("table");
            resetTab(tab);
            tab.style.border = "1px dashed #000000";
            toggleVisibility(this.className);
            tab.getElementsByTagName('th')[0].innerHTML = tab.ori_header +
                 "<br>Only showing reactions involving: " + this.innerHTML +
                 " (double-click to reset)";
        };
    }
};
var chempy_tabs = document.querySelectorAll('table.chempy_%(rsys_id)d');
[].forEach.call(chempy_tabs, function(tab){
    tab.ori_header = tab.getElementsByTagName('th')[0].innerHTML;
    tab.ondblclick = function(){
        resetTab(this);
        this.scrollIntoView();
    };
});
[].forEach.call(chempy_tabs, function(tab){
    resetTab(tab);
});
"""


def _js_rsys(cls_names_substances, substance_row_cls_irrel, rsys_id):
    # from https://stackoverflow.com/a/12786869/790973
    if cls_names_substances is None:
        return ''
    return _js_rsys_template % dict(
        cls_names_substances=json.dumps(cls_names_substances),
        substance_row_cls_irrel=json.dumps(substance_row_cls_irrel),
        rsys_id=rsys_id
    )


class JSPrinter(CSSPrinter):
    """ Prints javascript-enabled HTML representaions """

    def _print_ReactionSystem(self, rsys, **kwargs):
        tab = super(JSPrinter, self)._print_ReactionSystem(rsys, **kwargs)
        _script_tag = '<script type="text/javascript">%s</script>'
        substances = self._get('substances', **kwargs)
        cls_names_substances = list(map(_html_clsname, substances))
        return tab + _script_tag % _js_rsys(
            cls_names_substances=cls_names_substances,
            substance_row_cls_irrel={cns: [  # reactions not involving sk
                self._tr_id(rsys, i) for i in range(rsys.nr) if
                i not in rsys.substance_participation(sk)
            ] for cns, sk in zip(cls_names_substances, substances)},
            rsys_id=id(rsys)
        )


def javascript(obj, **settings):
    return JSPrinter(settings).doprint(obj)
