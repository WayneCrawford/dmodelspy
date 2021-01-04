#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test the functions/classes
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
# from future.builtins import *  # NOQA @UnusedWildImport

import inspect
import os
import unittest
import numpy as np
from dmodelspy import Sill
# import warnings
from itertools import cycle
import sys
sys.path.append('../')


class TestSillMethods(unittest.TestCase):
    """
    Test suite for dmodelspy operations.
    """
    def setUp(self):
        # self.path = os.path.dirname(os.path.abspath(inspect.getfile(
        #     inspect.currentframe())))
        # self.testing_path = os.path.join(self.path, "data")
        eps = 1e-8
        self.sill = Sill(0, 0, 1000, P_G=0.01, a=1000., nu=0.25)
        self.x = 1000. * np.array([eps, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75,
                              1., 1.5, 2., 3., 4., 5.])
        self.y = 2*self.x

    def test_sill_dV(self):
        """
        Test Sill calculated change in volume
        """
        self.assertAlmostEqual(self.sill.dV, 25723944.363136075)

    def test_sill_disp(self):
        """
        Test Sill calculated displacements
        """
        # Test for z=0
        z = 0
        u, v, w = self.sill.calc_displ(self.x, self.y, z)
        uref = [3.259e-08, 1.618e-01, 3.167e-01, 5.790e-01, 7.398e-01,
                7.746e-01, 6.990e-01, 3.627e-01, 1.593e-01, 3.800e-02,
                1.267e-02, 2.582e-03, 8.246e-04, 3.397e-04]
        vref = [6.519e-08, 3.236e-01, 6.334e-01, 1.158e+00, 1.480e+00,
                1.549e+00, 1.398e+00, 7.253e-01, 3.185e-01, 7.601e-02,
                2.534e-02, 5.165e-03, 1.649e-03, 6.795e-04]
        wref = [6.338e+00, 6.274e+00, 6.088e+00, 5.380e+00, 4.345e+00,
                3.189e+00, 2.138e+00, 6.295e-01, 1.875e-01, 2.722e-02,
                6.565e-03, 8.655e-04, 2.143e-04, 7.068e-05]
        for val, ref in zip(u, uref):
            self.assertAlmostEqual(val, ref, places=3)
        for val, ref in zip(v, vref):
            self.assertAlmostEqual(val, ref, places=3)
        for val, ref in zip(w, wref):
            self.assertAlmostEqual(val, ref, places=3)

        # Test for z=500
        z = 500
        u, v, w = self.sill.calc_displ(self.x, self.y, z)
        uref = [6.52369767e-09, 3.26048839e-02, 6.51687684e-02, 1.29955655e-01,
                1.85836745e-01, 1.97592812e-01, 1.42599260e-01, 4.76398821e-02,
                3.90532726e-02, 2.00031029e-02, 8.81152151e-03, 2.15762664e-03,
                6.89136874e-04, 3.61266220e-04]
        vref = [1.30473953e-08, 6.52097679e-02, 1.30337537e-01, 2.59911310e-01,
                3.71673491e-01, 3.95185625e-01, 2.85198521e-01, 9.52797643e-02,
                7.81065453e-02, 4.00062059e-02, 1.76230430e-02, 4.31525328e-03,
                1.37827375e-03, 7.22532440e-04]
        wref = [7.27274699e+00, 7.19989330e+00, 6.98381527e+00, 6.15198203e+00,
                4.87187786e+00, 3.35544954e+00, 2.02772724e+00, 5.42347042e-01,
                1.61948271e-01, 2.10103134e-02, 4.57399425e-03, 6.95581122e-05,
                6.58323848e-04, 2.98096757e-04]
        for val, ref in zip(u, uref):
            self.assertAlmostEqual(val, ref)
        for val, ref in zip(v, vref):
            self.assertAlmostEqual(val, ref)
        for val, ref in zip(w, wref):
            self.assertAlmostEqual(val, ref)

    def test_sill_tilt(self):
        """
        Test Sill calculated tilts
        """
        z = 0
        dwdx, dwdy = self.sill.calc_tilt(self.x, self.y, z)
        # print(dwdx, dwdy)
        dwdxref = [-1.01366737e-10, -5.03338383e-04, -9.85607834e-04, -1.80025941e-03,
                   -2.27030349e-03, -2.27559867e-03, -1.87457366e-03, -6.38664025e-04,
                   -1.69943950e-04, -1.77953776e-05, -3.26127806e-06, -2.70575019e-07,
                   -3.55031852e-08, -1.39514764e-08]
        dwdyref = [-2.02733475e-10, -1.00667675e-03, -1.97121557e-03, -3.60051851e-03,
                   -4.54060986e-03, -4.55121450e-03, -3.74918580e-03, -1.27736018e-03,
                   -3.39895936e-04, -3.55912253e-05, -6.52260969e-06, -5.41165109e-07,
                   -7.10161111e-08, -2.79001499e-08]
        for val, ref in zip(dwdx, dwdxref):
            self.assertAlmostEqual(val, ref)
        for val, ref in zip(dwdy, dwdyref):
            self.assertAlmostEqual(val, ref)

    def test_sill_strain(self):
        """
        Test Sill calculated tilts
        """
        z = 500
        eea, g1, g2 = self.sill.calc_strain(self.x, self.y, z)
        # print(eea, g1, g2)
        eea_ref = [ 1.30473826e-03,  1.30373737e-03,  1.30281874e-03,
                    1.28487886e-03,  1.03975915e-03,  2.44299413e-04,
                   -4.28086785e-04, -1.61933075e-05,  8.02265851e-06,
                   -1.85116933e-05, -8.28033274e-06, -6.20192219e-06,
                    4.07107591e-06,  2.35277399e-06]
        g1_ref = [ 0.00000000e+00,  2.74252965e-07,  3.33928467e-07,
                   8.80625980e-06,  1.19485084e-04,  4.46174228e-04,  
                   5.99073190e-04,  8.59432474e-05,  4.20501981e-05,  
                   2.71096902e-05,  1.02553796e-05,  4.58365417e-06,
                  -2.23541614e-06, -1.32461443e-06]
        g2_ref = [-2.03944689e-20, -3.65904259e-07, -4.45370174e-07,
                  -1.17361472e-05, -1.59297350e-04, -5.94909737e-04,
                  -7.98798413e-04, -1.14586535e-04, -5.60670125e-05,
                  -3.61459420e-05, -1.36733852e-05, -6.11270298e-06,
                   2.98154483e-06,  1.76688991e-06]
        for val, ref in zip(eea, eea_ref):
            self.assertAlmostEqual(val, ref)
        for val, ref in zip(g1, g1_ref):
            self.assertAlmostEqual(val, ref)
        for val, ref in zip(g2, g2_ref):
            self.assertAlmostEqual(val, ref)


def suite():
    return unittest.makeSuite(TestSillMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
