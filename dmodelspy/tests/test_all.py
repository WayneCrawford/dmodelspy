#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test the functions/classes
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
# from future.builtins import *  # NOQA @UnusedWildImport

# import inspect
# import os
import dataclasses
import unittest
import numpy as np
from dmodelspy import Sill, Sphere3D, Spheroid
# import warnings
# from itertools import cycle
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
        self.model = Sill(0, 0, 1000, P_G=0.01, a=1000., nu=0.25)
        self.x = 1000. * np.array([eps, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75,
                                   1., 1.5, 2., 3., 4., 5.])
        self.y = 2*self.x

    def test_dV(self):
        """
        Test Sill calculated change in volume
        """
        self.assertAlmostEqual(self.model.dV, 25723944.363136075)

    def test_displ(self):
        """
        Test Sill calculated displacements
        """
        # Test for z=0
        z = 0
        u, v, w = self.model.calc_displ(self.x, self.y, z)
        # print(f'{u=}\n{v=}\n{w=}')
        uref = [
            3.25948703e-08, 1.61813815e-01, 3.16694648e-01, 5.78977374e-01,
            7.39846808e-01, 7.74632048e-01, 6.99022455e-01, 3.62654888e-01,
            1.59265750e-01, 3.80027255e-02, 1.26678682e-02, 2.58232749e-03,
            8.24597350e-04, 3.39725418e-04]
        vref = [
            6.51897405e-08, 3.23627631e-01, 6.33389296e-01, 1.15795475e+00,
            1.47969362e+00, 1.54926410e+00, 1.39804491e+00, 7.25309775e-01,
            3.18531499e-01, 7.60054511e-02, 2.53357364e-02, 5.16465498e-03,
            1.64919470e-03, 6.79450837e-04]
        wref = [
            6.33750037e+00, 6.27436367e+00, 6.08758210e+00, 5.38031833e+00,
            4.34494951e+00, 3.18850224e+00, 2.13827971e+00, 6.29510464e-01,
            1.87537770e-01, 2.72191189e-02, 6.56480277e-03, 8.65511170e-04,
            2.14267298e-04, 7.06838146e-05]
        for val, ref in zip(u, uref):
            self.assertAlmostEqual(val, ref)
        for val, ref in zip(v, vref):
            self.assertAlmostEqual(val, ref)
        for val, ref in zip(w, wref):
            self.assertAlmostEqual(val, ref)

        # Test for z=500
        z = 500
        u, v, w = self.model.calc_displ(self.x, self.y, z)
        # print(f'{u=}\n{v=}\n{w=}')
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

    def test_tilt(self):
        """
        Test calculated tilts
        """
        z = 0
        dwdx, dwdy = self.model.calc_tilt(self.x, self.y, z)
        # print(dwdx, dwdy)
        dwdxref = [
            -1.01366737e-10, -5.03338383e-04, -9.85607834e-04, -1.80025941e-03,
            -2.27030349e-03, -2.27559867e-03, -1.87457366e-03, -6.38664025e-04,
            -1.69943950e-04, -1.77953776e-05, -3.26127806e-06, -2.70575019e-07,
            -3.55031852e-08, -1.39514764e-08]
        dwdyref = [
            -2.02733475e-10, -1.00667675e-03, -1.97121557e-03, -3.60051851e-03,
            -4.54060986e-03, -4.55121450e-03, -3.74918580e-03, -1.27736018e-03,
            -3.39895936e-04, -3.55912253e-05, -6.52260969e-06, -5.41165109e-07,
            -7.10161111e-08, -2.79001499e-08]
        for val, ref in zip(dwdx, dwdxref):
            self.assertAlmostEqual(val, ref)
        for val, ref in zip(dwdy, dwdyref):
            self.assertAlmostEqual(val, ref)

    def test_strain(self):
        """
        Test calculated tilts
        """
        z = 500
        eea, g1, g2 = self.model.calc_strain(self.x, self.y, z)
        # print(eea, g1, g2)
        eea_ref = [1.30473826e-03,  1.30373737e-03,  1.30281874e-03,
                   1.28487886e-03,  1.03975915e-03,  2.44299413e-04,
                   -4.28086785e-04, -1.61933075e-05,  8.02265851e-06,
                   -1.85116933e-05, -8.28033274e-06, -6.20192219e-06,
                   4.07107591e-06,  2.35277399e-06]
        g1_ref = [0.00000000e+00,  2.74252965e-07,  3.33928467e-07,
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


class TestSpheroidMethods(unittest.TestCase):
    """
    Test suite for dmodelspy operations.
    """
    def setUp(self):
        # self.path = os.path.dirname(os.path.abspath(inspect.getfile(
        #     inspect.currentframe())))
        # self.testing_path = os.path.join(self.path, "data")
        eps = 1e-8
        self.model = Spheroid(0, 0, 1000, P_G=0.1, a=1000., nu=0.25,
                              asrat=0.5, mu=9.6e9, phi=30, theta=90)
        self.x = 1000. * np.array([eps, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75,
                                   1., 1.5, 2., 3., 4., 5.])
        self.y = 2*self.x

    def test_dV(self):
        """
        Test Sill calculated change in volume
        """
        # print(f'{self.model.dV=}')
        self.assertAlmostEqual(self.model.dV, 78539816.33974484)

    def test_displ(self):
        """
        Test Sill calculated displacements
        """
        # Test for z=0
        z = 0
        u, v, w = self.model.calc_displ(self.x, self.y, z)
        # print(f'{u=}\n{v=}\n{w=}')
        uref = [0.00611681, 3.79838517, 4.92779013, 4.86514901, 4.26886095,
                3.65548878, 3.11046954, 2.09104458, 1.45043479, 0.78000386,
                0.4747821 , 0.22441055, 0.12914373, 0.08354751]
        vref = [0.0105949 , 7.59556548, 9.85487406, 9.73002142, 8.53759228,
                7.31091104, 6.22090155, 4.18207816, 2.90086334, 1.56000404,
                0.94956306, 0.44882536, 0.25827517, 0.16708555]
        wref = [
            2.83795239e+01, 2.43460713e+01, 1.90438377e+01, 1.24551547e+01,
            8.69674299e+00, 6.28370297e+00, 4.64379395e+00, 2.34870503e+00,
            1.29893265e+00, 4.92093842e-01, 2.29791142e-01, 7.36910005e-02,
            3.20088105e-02, 1.66160410e-02]
        for val, ref in zip(u, uref):
            self.assertAlmostEqual(val, ref)
        for val, ref in zip(v, vref):
            self.assertAlmostEqual(val, ref)
        for val, ref in zip(w, wref):
            self.assertAlmostEqual(val, ref)

        # Test for asrat=0.75
        newmodel = dataclasses.replace(self.model)
        newmodel.asrat = 0.75
        u, v, w = newmodel.calc_displ(self.x, self.y, z)
        # print(f'{u=}\n{v=}\n{w=}')
        uref = [
            3.49799101e-03, 4.04026383e+00, 6.88839922e+00, 8.94360197e+00,
            8.63655430e+00, 7.65899931e+00, 6.59106535e+00, 4.41203510e+00,
            3.02081849e+00, 1.59351883e+00, 9.60050644e-01, 4.49685776e-01,
            2.57862421e-01, 1.66501357e-01]
        vref = [
            6.05894216e-03, 8.07977642e+00, 1.37761960e+01, 1.78868826e+01,
            1.72729433e+01, 1.53179114e+01, 1.31820856e+01, 8.82405890e+00,
            6.04164028e+00, 3.18703082e+00, 1.92008153e+00, 8.99400156e-01,
            5.15722611e-01, 3.33036027e-01]
        wref = [
            5.09365156e+01, 4.86176669e+01, 4.31055658e+01, 3.08682177e+01,
            2.16382944e+01, 1.53481771e+01, 1.10796031e+01, 5.31325343e+00,
            2.83034909e+00, 1.02859085e+00, 4.71004463e-01, 1.48600099e-01,
            6.41435611e-02, 3.31974397e-02]
        for val, ref in zip(u, uref):
            self.assertAlmostEqual(val, ref)
        for val, ref in zip(v, vref):
            self.assertAlmostEqual(val, ref)
        for val, ref in zip(w, wref):
            self.assertAlmostEqual(val, ref)

    def test_tilt(self):
        """
        Test calculated tilts
        """
        z = 0
        dwdx, dwdy = self.model.calc_tilt(self.x, self.y, z)
        # print(f'{dwdx=}, {dwdy=}')
        dwdxref = [
            -5.31247591e-05, -2.34633563e-02, -1.80749962e-02, -9.56594029e-03,
            -5.89842024e-03, -3.92833769e-03, -2.72083428e-03, -1.19214994e-03,
            -5.73990835e-04, -1.67805032e-04, -6.26358019e-05, -1.40920871e-05,
            -4.67950450e-06, -1.96108158e-06]
        dwdyref = [
            -9.20172214e-05, -4.69380483e-02, -3.61583016e-02, -1.91337189e-02,
            -1.17975069e-02, -7.85699542e-03, -5.44184515e-03, -2.38435235e-03,
            -1.14800062e-03, -3.35613625e-04, -1.25272504e-04, -2.81842237e-05,
            -9.35900900e-06, -3.92214134e-06]
        for val, ref in zip(dwdx, dwdxref):
            self.assertAlmostEqual(val, ref)
        for val, ref in zip(dwdy, dwdyref):
            self.assertAlmostEqual(val, ref)

    def test_strain(self):
        """
        Test calculated tilts
        """
        z = 500
        eea, g1, g2 = self.model.calc_strain(self.x, self.y, z)
        # print(f'{eea=}, {g1=}, {g2=}')
        eea_ref = [
            6.576667202e+02, -1.76250313e-02, -9.43997550e-03, -1.98817454e-03,
            -2.35329097e-04, -4.73239338e-05, -1.79024227e-04, -4.31972049e-04,
            -4.36852864e-04, -2.80061961e-04, -1.62380354e-04, -6.24405064e-05,
            -2.93822946e-05, -1.59880868e-05]
        g1_ref = [
            -1.00230098e-01,  7.84566251e-01,  1.93398686e-01,  4.60957796e-02,
            1.98264460e-02,  1.10800890e-02,  7.15750151e-03,  3.27905396e-03,
            1.84347925e-03,  7.52603976e-04,  3.70551218e-04,  1.24790845e-04,
            5.68674798e-05,  3.03632524e-05]
        g2_ref = [
            -1.73597226e-01, -1.04629529e+00, -2.57888105e-01, -6.14630927e-02,
            -2.64358131e-02, -1.47737771e-02, -9.54304827e-03, -4.37203063e-03,
            -2.45808090e-03, -1.00354482e-03, -4.93809922e-04, -1.67495253e-04,
            -7.41415292e-05, -3.95819271e-05]
        for val, ref in zip(eea, eea_ref):
            self.assertAlmostEqual(val, ref)
        for val, ref in zip(g1, g1_ref):
            self.assertAlmostEqual(val, ref)
        for val, ref in zip(g2, g2_ref):
            self.assertAlmostEqual(val, ref)


class TestSphere3DMethods(unittest.TestCase):
    """
    Test suite for Sphere3D
    """
    def setUp(self):
        eps = 1e-8
        self.model = Sphere3D(0, 0, 1000, P_G=0.001273, a=500., nu=0.25)
        self.x = 1000. * np.array([eps, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75,
                                   1., 1.5, 2., 3., 4., 5.])
        self.y = 2*self.x

    def test_dV(self):
        """
        Test Sill calculated change in volume
        """
        # print(f'{self.model.dV=}')
        self.assertAlmostEqual(self.model.dV, 499905.9310024759)

    def test_displ(self):
        """
        Test Sill calculated displacements
        """
        # Test for z=0
        z = 0
        u, v, w = self.model.calc_displ(self.x, self.y, z)
        # print(f'{u=}\n{v=}\n{w=}')
        uref = [
            1.32943484e-09, 6.51604186e-03, 1.22944738e-02, 1.98697001e-02,
            2.20821255e-02, 2.09589081e-02, 1.84978860e-02, 1.22964007e-02,
            8.20167945e-03, 4.16954478e-03, 2.46545811e-03, 1.13669114e-03,
            6.47782407e-04, 4.17106439e-04]
        vref = [
            2.65886969e-09, 1.30320837e-02, 2.45889475e-02, 3.97394002e-02,
            4.41642509e-02, 4.19178163e-02, 3.69957721e-02, 2.45928015e-02,
            1.64033589e-02, 8.33908956e-03, 4.93091622e-03, 2.27338228e-03,
            1.29556481e-03, 8.34212879e-04]
        wref = [
            1.38345543e-01, 1.35570150e-01, 1.27771930e-01, 1.02903277e-01,
            7.59233491e-02, 5.38290043e-02, 3.78762014e-02, 1.66873135e-02,
            8.32116994e-03, 2.81197902e-03, 1.24556546e-03, 3.82493579e-04,
            1.63428422e-04, 8.41719449e-05]
        for val, ref in zip(u, uref):
            self.assertAlmostEqual(val, ref, places=3)
        for val, ref in zip(v, vref):
            self.assertAlmostEqual(val, ref, places=3)
        for val, ref in zip(w, wref):
            self.assertAlmostEqual(val, ref, places=3)

        # Test for z=500
        z = 500
        u, v, w = self.model.calc_displ(self.x, self.y, z)
        # print(f'{u=}\n{v=}\n{w=}')
        uref = [
            3.20137100e-09, 1.49321871e-02, 2.46228584e-02, 2.71014316e-02,
            2.13904775e-02, 1.63299312e-02, 1.29460289e-02, 8.46077472e-03,
            6.10944137e-03, 3.52022850e-03, 2.22063876e-03, 1.08150656e-03,
            6.29468938e-04, 4.09435896e-04]
        vref = [
            6.40274201e-09, 2.98643742e-02, 4.92457167e-02, 5.42028631e-02,
            4.27809550e-02, 3.26598625e-02, 2.58920577e-02, 1.69215494e-02,
            1.22188827e-02, 7.04045700e-03, 4.44127752e-03, 2.16301311e-03,
            1.25893788e-03, 8.18871792e-04]
        wref = [
            2.51540128e-01, 2.40993373e-01, 2.14246492e-01, 1.48406552e-01,
            9.59289851e-02, 6.14467410e-02, 4.01094832e-02, 1.58192585e-02,
            7.48059876e-03, 2.41767141e-03, 1.05265791e-03, 3.19131803e-04,
            1.35728852e-04, 6.97552126e-05]
        for val, ref in zip(u, uref):
            self.assertAlmostEqual(val, ref)
        for val, ref in zip(v, vref):
            self.assertAlmostEqual(val, ref)
        for val, ref in zip(w, wref):
            self.assertAlmostEqual(val, ref)

    def test_tilt(self):
        """
        Test Sill calculated tilts
        """
        z = 0
        dwdx, dwdy = self.model.calc_tilt(self.x, self.y, z)
        # print(f'{dwdx=}\n{dwdy=}')
        dwdxref = [
            -4.51511329e-12, -2.18334823e-05, -3.96006579e-05, -5.53909492e-05,
            -5.02287422e-05, -3.78609768e-05, -2.63971223e-05, -1.01258966e-05,
            -4.23716855e-06, -1.04251473e-06, -3.57812582e-07, -7.50230295e-08,
            -2.42461139e-08, -1.00296863e-08]
        dwdyref = [
            -9.03022936e-12, -4.36670180e-05, -7.92016755e-05, -1.10783424e-04,
            -1.00459586e-04, -7.57237533e-05, -5.27954811e-05, -2.02521546e-05,
            -8.47444376e-06, -2.08504339e-06, -7.15628035e-07, -1.50046340e-07,
            -4.84922798e-08, -2.00593864e-08]
        for val, ref in zip(dwdx, dwdxref):
            self.assertAlmostEqual(val, ref)
        for val, ref in zip(dwdy, dwdyref):
            self.assertAlmostEqual(val, ref)

    def test_strain(self):
        """
        Test Sill calculated tilts
        """
        z = 500
        eea, g1, g2 = self.model.calc_strain(self.x, self.y, z)
        # print(f'{eea=}\n{g1=}\n{g2=}')
        eea_ref = [
            6.40183445e-04,  5.56525199e-04,  3.72814830e-04,  9.33357298e-05,
            1.23708367e-05, -7.58975525e-07, -1.28058277e-06, -7.90456525e-07,
            -1.29776909e-06, -1.19351865e-06, -7.61512743e-07, -3.05192563e-07,
            -1.43320722e-07, -7.71327085e-08]
        g1_ref = [
            1.08420217e-18, 2.44086919e-05, 7.17510670e-05, 1.06598641e-04,
            7.81385382e-05, 4.94453716e-05, 3.18390034e-05, 1.40115472e-05,
            8.10999971e-06, 3.53229677e-06, 1.78929205e-06, 6.15718394e-07,
            2.74833180e-07, 1.44544264e-07]
        g2_ref = [
            -1.45192602e-18, -3.25515964e-05, -9.56836280e-05, -1.42141285e-04,
            -1.04186996e-04, -6.59275083e-05, -4.24520148e-05, -1.86820598e-05,
            -1.08133362e-05, -4.70972902e-06, -2.38572238e-06, -8.20957742e-07,
            -3.66444203e-07, -1.92725670e-07]
        for val, ref in zip(eea, eea_ref):
            self.assertAlmostEqual(val, ref)
        for val, ref in zip(g1, g1_ref):
            self.assertAlmostEqual(val, ref)
        for val, ref in zip(g2, g2_ref):
            self.assertAlmostEqual(val, ref)


def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSpheroidMethods, 'test'))
    suite.addTest(unittest.makeSuite(TestSillMethods, 'test'))
    suite.addTest(unittest.makeSuite(TestSphere3DMethods, 'test'))
    return suite
    # return unittest.makeSuite(TestSphere3DMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
