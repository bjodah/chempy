from chempy import Reaction
from chempy.util.testing import requires


@requires('sympy')
def test_Reaction_string():
    from sympy import S
    r = Reaction({'A': 1, 'B': 2}, {'C': S(3)/2}, checks=[
        chk for chk in Reaction.default_checks if chk != 'all_integral'])
    assert r.string() == 'A + 2 B -> 3/2 C'
