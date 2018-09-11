from nose.tools import *
import numpy as np

from nanosims_analysis.data_structures import IsotopeData

class TestClass:

    @classmethod
    def setup_class(cls):
        test_data = np.array([[[0.92672997, 0.98710138, 0.69649032],
                               [0.13295207, 0.66058395, 0.71824247],
                               [0.39265736, 0.03355262, 0.43079679]],

                              [[0.0554341 , 0.59301986, 0.43093438],
                               [0.46924799, 0.50525259, 0.38235398],
                               [0.85084851, 0.63285618, 0.93021085]]])
        cls.testIsotope = IsotopeData("test", test_data)

    def test_dummy(self):
        dwell_time = 0.006
        dead_time  = 4*np.power(10.0,-9)
        corrected = np.array([[[154.45509043, 164.51700493, 116.0817739 ],
                               [ 22.1586803 , 110.09737349, 119.70713565],
                               [ 65.44291046,   5.59210346,  71.79948562]],

                              [[  9.23901701,  98.83668241,  71.8224173 ],
                               [ 78.2080228 ,  84.20879336,  63.72567958],
                               [141.80816544, 105.4760745 , 155.03523781]]])
        self.testIsotope.perform_deadtime_correction(dwell_time, dead_time)
        assert_true(np.allclose(self.testIsotope._data,corrected))
                                
