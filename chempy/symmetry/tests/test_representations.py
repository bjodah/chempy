# -*- coding: utf-8 -*-

import pytest
import numpy as np
from ..representations import Reducible

# Test data with (reducible representation, num of irreducibles, point group)
reducible_data = (
            ((2, 0), (1, 1), 'cs'),
            ((3, 1), (2, 1), 'ci'),
            ((2), (2), 'c1'),
            ((3, 1), (2, 1), 'c2'),
            ((3, 0, 0), (1, 1), 'C3'),
            ((4, 0, 0, 0), (1, 1, 1), 'c4'),
            ((21, 1, 1, 1, 1), (5, 4, 4), 'c5'),
            ((6, 0, 0, 0, 0, 0), (1, 1, 1, 1), 'c6'),
            ([9, -1, 3, 1], [3, 1, 3, 2], 'c2v'),
            ([3, 0, 1], (1, 0, 1), 'c3v'),
            ((6, 0, 2, 0, 0), (1, 1, 1, 1, 1), 'c4v'),
            ((4, 0, 2, 2), (2, 1, 0, 1), 'c2h'),
            ((15, 0, 0, 7, -2, -2), (3, 4, 2, 1), 'c3h'),
            ((5, -1, 1, -1, -1, 1, -5, 1), (0, 0, 1, 1, 2, 0), 'c4h'),
            ((4, 0, 0, 0), (1, 1, 1, 1), 'd2'),
            ((5, -1, 1, 3, -1), (1, 0, 2, 0, 1), 'd2d'),
            ((5, 2, 1, 3, 0, 3), (2, 0, 1, 0, 1, 0), 'd3h'),
            ((4, 0, 0, 2, 0, 0, 0, 4, 2, 0), (1, 0, 1, 0, 0, 0, 0, 0, 0, 1), 'd4h'),
            ((6, 0, 0, 0, -2, 0, 0, 0, 0, -6, 0, 2), (0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1), 'd6h'),
            ([4, 1, 0, 0, 2], (1, 0, 0, 0, 1), 'Td'),
            ([8, -1, 4, 0, -2], (0, 1, 2, 1, 0), 'td'),
            ((6, 0, 0, 2, 2, 0, 0, 0, 4, 2), (1, 0, 1, 0, 0, 0, 0, 0, 1, 0), 'Oh'),
            ((7, 1, 1, 3, 3, 1, 1, 1, 5, 3), (2, 0, 1, 0, 0, 0, 0, 0, 1, 0), 'oh'),
            )


@pytest.mark.parametrize('gamma, true_n_irred, group', reducible_data)
def test_decomp(gamma, true_n_irred, group):
    calc_n_irred = Reducible(gamma, group).decomp()

    assert np.all(true_n_irred == calc_n_irred)

all_motion_data = (
    ((6, 0, 0, 6), (18, 0, 0, 6), 'c2h'),
    ((3, 1, 3, 1), (9, -1, 3, 1), 'c2v'),
    ((4, 1, 2), (12, 0, 2), 'c3v'),
    ((6, 2, 2, 4, 2), (18, 2, -2, 4, 2), 'c4v'),
    ((8, 2, 0, 0, 0, 4), (24, 0, 0, 0, 0, 4), 'd3d'),
    ((6, 3, 2, 4, 1, 4), (18, 0, -2, 4, -2, 4), 'd3h'),
    ((5, 1, 1, 3, 1, 1, 1, 5, 3, 1), (15, 1, -1, -3, -1, -3, -1, 5, 3, 1), 'd4h'),
    ((5, 2, 1, 1, 3), (15, 0, -1, -1, 3), 'td'),
    ((7, 1, 1, 3, 3, 1, 1, 1, 5, 3), (21, 0, -1, 3, -3, -3, -1, 0, 5, 3), 'oh')
    )

@pytest.mark.parametrize('stationary_atoms, reducible_rep, group', all_motion_data)
def test_from_atoms(stationary_atoms, reducible_rep, group):

    assert np.all(Reducible.from_atoms(stationary_atoms, group).gamma ==
                  np.array(reducible_rep))


class Test_ReducibleMethods():

    def test_decomp(self):
        water = Reducible([9, -1, 3, 1], 'c2v', all_motion=True)
        # trans-dichloroethene
        tDCE = Reducible([18, 0, 0, 6], 'C2h', all_motion=True)
        # cis-dichloroethene
        cDCE = Reducible([18, 0, 6, 0], 'C2v', all_motion=True)
        # bromopentacarbonylmolybdenum(0)
        pCOBrMn = Reducible([5, 1, 1, 3, 1], 'C4v', all_motion=False)

        assert np.all(water.decomp() == np.array([3, 1, 3, 2]))
        assert np.all(tDCE.decomp() == np.array([6, 3, 3, 6]))
        assert np.all(cDCE.decomp() == np.array([6, 3, 6, 3]))
        assert np.all(pCOBrMn.decomp() == np.array([2, 0, 1, 0, 1]))
        assert np.all(pCOBrMn.decomp(to_dict=True) ==
                      {'A1': 2, 'A2': 0, 'B1': 1, 'B2': 0, 'E': 1})

    def test_vibe(self):
        water = Reducible([9, -1, 3, 1], 'c2v', all_motion=True)
        tDCE = Reducible([18, 0, 0, 6], 'C2h', all_motion=True)

        assert np.all(water.vibe_modes() == np.array([2, 0, 1, 0]))
        assert np.all(tDCE.vibe_modes() == np.array([5, 1, 2, 4]))
        assert np.all(tDCE.vibe_modes(to_dict=True) ==
                      {'Ag': 5, 'Bg': 1, 'Au': 2, 'Bu': 4})

    def test_ir(self):
        water = Reducible([9, -1, 3, 1], 'c2v', all_motion=True)
        pentCOMn = Reducible([5, 2, 1, 3, 0, 3], 'd3h', all_motion=False)
        tDCE = Reducible([18, 0, 0, 6], 'c2h', all_motion=True)
        pCOBrMn = Reducible([5, 1, 1, 3, 1], 'C4v', all_motion=False)

        assert np.all(water.ir_active() == np.array([2, 0, 1, 0]))
        assert np.all(pentCOMn.ir_active() == np.array([0, 0, 1, 0, 1, 0]))
        assert np.all(tDCE.ir_active() == np.array([0, 0, 2, 4]))
        assert np.all(pCOBrMn.ir_active() == np.array([2, 0, 0, 0, 1]))
        assert np.all(water.ir_active(to_dict=True) ==
                      {'A1': 2, 'A2': 0, 'B1': 1, 'B2': 0})

    def test_raman(self):
        water = Reducible([9, -1, 3, 1], 'c2v', all_motion=True)
        tDCE = Reducible([18, 0, 0, 6], 'c2h', all_motion=True)
        pCOBrMn = Reducible([5, 1, 1, 3, 1], 'C4v', all_motion=False)

        assert np.all(water.raman_active() == np.array([2, 0, 1, 0]))
        assert np.all(tDCE.raman_active() == np.array([5, 1, 0, 0]))
        assert np.all(pCOBrMn.raman_active() == np.array([2, 0, 1, 0, 1]))
        assert np.all(tDCE.raman_active(to_dict=True) ==
                      {'Ag': 5, 'Bg': 1, 'Au': 0, 'Bu': 0})


    def test_from_irred(self):
        test_rep = Reducible.from_irred([1, 0, 1, 0], 'c2v')
        true_rep = Reducible([2, 0, 2, 0], 'c2v', all_motion=False)
        assert np.all(test_rep.gamma == true_rep.gamma)

    @pytest.mark.xfail()
    def test_fail(self):
        water = Reducible([9, -1, 3, 1], 'c2v', all_motion=False)
        assert np.all(water.ir_active() == np.array([2, 0, 1, 0]))
