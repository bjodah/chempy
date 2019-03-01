import math
from ..rendering import eval_template

def test_eval_template():
    assert eval_template("2 OH -> H2O2; ${6*pi*arg}/M/s", arg=1/math.pi) == "2 OH -> H2O2; 6.0/M/s"
