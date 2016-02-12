import numpy as np
from chempy.tests.ammonical_cupric_solution import get_ammonical_cupric_eqsys


class TimeEqsys:

    def setup(self):
        self.eqsys, self.c0 = get_ammonical_cupric_eqsys()

    def time_roots(self):
        x, new_inits, success = self.eqsys.roots(
            self.c0, np.logspace(-3, 0, 50), 'NH3')
        assert all(success)

    def time_roots_symengine(self):
        from symengine import Lambdify
        x, new_inits, success = self.eqsys.roots(
            self.c0, np.logspace(-3, 0, 50), 'NH3',
            lambdify=Lambdify, lambdify_unpack=False)
        assert all(success)

    def time_roots_no_propagate(self):
        x, new_inits, success = self.eqsys.roots(
            self.c0, np.logspace(-3, 0, 50), 'NH3', propagate=False)
        assert all(success)


if __name__ == '__main__':
    import time
    te = TimeEqsys()
    te.setup()

    # t1 = time.time()
    # te.time_roots_symengine()
    # print(time.time()-t1)

    t1 = time.time()
    te.time_roots()
    print(time.time()-t1)
