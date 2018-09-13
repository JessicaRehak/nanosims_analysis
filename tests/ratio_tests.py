from nose.tools import *
import numpy as np

from nanosims_analysis.importer import Importer
from nanosims_analysis.data_structures import IsotopeData
from nanosims_analysis.data_structures import RatioData

class TestClass:

    @classmethod
    def setup_class(cls):
        cls.numerator_isotope = IsotopeData("numerator",
                                             np.array([[[782, 645, 275],
                                                        [733, 875, 590],
                                                        [902, 708, 739]],

                                                       [[789, 433, 716],
                                                        [ 72, 367, 390],
                                                        [465,  58, 995]]]))


        cls.denominator_isotope = IsotopeData("denominator",
                                              np.array([[[500, 475, 370],
                                                         [779,  37, 769],
                                                         [ 33, 796, 741]],
                                            
                                                        [[502, 618, 573],
                                                         [963, 473, 858],
                                                         [ 31, 734,   0]]]))

    def test_ratio_init(self):
        testRatio = RatioData("test_ratio",
                              self.numerator_isotope,
                              self.denominator_isotope)

        ans = np.array([[[ 1.564     ,  1.35789474,  0.74324324],
                         [ 0.94094994, 23.64864865,  0.76723017],
                         [27.33333333,  0.88944724,  0.99730094]],

                        [[ 1.57171315,  0.70064725,  1.2495637 ],
                         [ 0.07476636,  0.77589852,  0.45454545],
                         [15.        ,  0.07901907,         0]]])
        assert_true(np.allclose(ans, testRatio.get_data()))

    @raises(RuntimeError)
    def test_ratio_deadtime(self):
        testRatio = RatioData("test_ratio",
                              self.numerator_isotope,
                              self.denominator_isotope)
        testRatio.perform_deadtime_correction(2.0, 3.0)
