from nose.tools import *
import numpy as np

from nanosims_analysis.data_structures import IsotopeData

class TestClass:

    @classmethod
    def setup_class(cls):
        test_data = [[[0.92672997, 0.98710138, 0.69649032],
                      [0.13295207, 0.66058395, 0.71824247],
                      [0.39265736, 0.03355262, 0.43079679]],

                     [[0.0554341 , 0.59301986, 0.43093438],
                      [0.46924799, 0.50525259, 0.38235398],
                      [0.85084851, 0.63285618, 0.93021085]]]
        cls.testIsotope = IsotopeData("test", test_data)

    def test_dummy(self):
        # dwell_time = 0.006
        # dead_time  = 4*10**-9
        assert_true(True)
