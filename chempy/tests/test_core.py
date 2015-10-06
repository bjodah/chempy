from ..core import ActivityProduct, ionic_strength


def test_ionic_strength():
    assert abs(ionic_strength([0.1, 1.3, 2.1, 0.7],
                              [-1, 2, -3, 4]) - 17.7) < 1e-14


def test_ActivityProduct():
    ap = ActivityProduct((1, -2, 3), 17, 42)
    assert ap.stoich == (1, -2, 3)
    assert ap.args == (17, 42)
