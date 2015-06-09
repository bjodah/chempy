from ..henry import Henry


def test_Henry():
    kH_O2 = Henry(1.2e-3, 1800, ref='carpenter_1966')
    assert abs(kH_O2.get_kH_at_T(298.15) - 1.2e-3) < 1e-4
    assert abs(kH_O2.get_c_at_T_and_P(290, 1) - 0.001421892) < 1e-8
    assert abs(kH_O2.get_P_at_T_and_c(310, 1e-3) - 1.05) < 1e-3

    try:
        import quantities as pq
        import numpy as np
        from .util import allclose

        dm3 = 1e-3*pq.metre**3
        molar = pq.mol / dm3
        kH_H2 = Henry(7.8e-4*molar/pq.atm, 640*pq.K, units=pq, ref='dean_1992')

        assert allclose(kH_H2.get_kH_at_T(
            300*pq.K), 7.697430323e-4*molar/pq.atm)
        kH = kH_H2.get_c_at_T_and_P(
            np.linspace(297.5, 298.65, 3)*pq.K, .1*pq.bar)
        assert allclose(kH, 7.7e-5*molar, rtol=1e-5, atol=1e-6*molar)
        kH = kH_H2.get_P_at_T_and_c(
            298.15*pq.K, np.linspace(2e-3, 2.1e-3, 3)*molar)
        assert allclose(kH, 2.65*pq.atm, rtol=1e-5, atol=0.2*pq.bar)
    except ImportError:
        pass
