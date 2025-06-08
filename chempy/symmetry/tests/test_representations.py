# -*- coding: utf-8 -*-

import pytest
import numpy as np
from ..representations import Reducible

# Test data with (reducible representation, num of irreducibles, point group)
test_data = (
            ((2), (2), 'c1'),
            ((2, 0), (1, 1), 'cs'),
            ((3, 1), (2, 1), 'ci'),
            ((3, 1), (2, 1), 'c2'),
            ((3, 0, 0), (1, 1), 'C3'),
            ((4, 0, 0, 0), (1, 1, 1), 'c4'),
            ((21, 1, 1, 1, 1), (5, 4, 4), 'c5'),
            ([9, -1, 3, 1], [3, 1, 3, 2], 'c2v'),
            ([3, 0, 1], (1, 0, 1), 'c3v'),
            ((6, 0, 2, 0, 0), (1, 1, 1, 1, 1), 'c4v'),
            ((4, 0, 2, 2), (2, 1, 0, 1), 'c2h'),
            ([15, 0, 0, 7, -2, -2], [3, 4, 2, 1], 'c3h'),
            ((4, 0, 0, 0), (1, 1, 1, 1), 'd2'),
            ([4, 1, 0, 0, 2], (1, 0, 0, 0, 1), 'Td'),
            ([4, 0, 0, 2, 0, 0, 0, 4, 2, 0],
             (1, 0, 1, 0, 0, 0, 0, 0, 0, 1), 'd4h'),
            ((5, 2, 1, 3, 0, 3), (2, 0, 1, 0, 1, 0), 'd3h'),
            ((6, 0, 0, 0, -2, 0, 0, 0, 0, -6, 0, 2),
             (0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1), 'd6h'),
            ([6, 0, 0, 2, 2, 0, 0, 0, 4, 2],
             (1, 0, 1, 0, 0, 0, 0, 0, 1, 0), 'Oh'),
            ((6, 0, 0, 0, 0, 0), (1, 1, 1, 1), 'c6')
            )


@pytest.mark.parametrize('gamma, true_n_irred, group', test_data)
def test_decomp(gamma, true_n_irred, group):
    calc_n_irred = Reducible(gamma, group).decomp()

    assert np.all(true_n_irred == calc_n_irred)


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


    def test_from_atoms(self):
        assert np.all(Reducible.from_atoms([6, 0, 0, 6], 'c2h').gamma ==
                      np.array([18, 0, 0, 6]))

        assert np.all(Reducible.from_atoms([3, 1, 3, 1], 'c2v').gamma ==
                      np.array([9, -1, 3, 1]))

        assert np.all(Reducible.from_atoms([21, 1, 1, 1, 1, 1, 1, 5], 'd5h').
                      gamma == np.array([63, 3, -1, -1, 1, 1, -2, 5]))

    def test_from_irred(self):
        test_rep = Reducible.from_irred([1, 0, 1, 0], 'c2v')
        true_rep = Reducible([2, 0, 2, 0], 'c2v', all_motion=False)
        assert np.all(test_rep.gamma == true_rep.gamma)

    @pytest.mark.xfail()
    def test_fail(self):
        water = Reducible([9, -1, 3, 1], 'c2v', all_motion=False)
        assert np.all(water.ir_active() == np.array([2, 0, 1, 0]))
