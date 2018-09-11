from nose.tools import *
import numpy as np

from nanosims_analysis.importer import Importer
from nanosims_analysis.data_structures import IsotopeData

class TestClass:

    @classmethod
    def setup_class(cls):
        cls.test_data = np.array([[[0.92672997, 0.98710138, 0.69649032],
                               [0.13295207, 0.66058395, 0.71824247],
                               [0.39265736, 0.03355262, 0.43079679]],

                              [[0.0554341 , 0.59301986, 0.43093438],
                               [0.46924799, 0.50525259, 0.38235398],
                               [0.85084851, 0.63285618, 0.93021085]]])
        cls.dwell_time = 0.006
        cls.dead_time  = 4*np.power(10.0,-9)
        cls.dt_corrected = np.array([[[154.45509043, 164.51700493, 116.0817739 ],
                                      [ 22.1586803 , 110.09737349, 119.70713565],
                                      [ 65.44291046,   5.59210346,  71.79948562]],

                                     [[  9.23901701,  98.83668241,  71.8224173 ],
                                      [ 78.2080228 ,  84.20879336,  63.72567958],
                                      [141.80816544, 105.4760745 , 155.03523781]]])

    def test_isotope_data_deadtime(self):
        testIsotope = IsotopeData("test", self.test_data)
        testIsotope.perform_deadtime_correction(dwell_time = self.dwell_time,
                                                dead_time = self.dead_time)
        #print(self.testIsotope._data)
        assert_true(np.allclose(testIsotope._data,self.dt_corrected))

    @raises(RuntimeError)
    def test_isotope_data_double_deadtime(self):
        testIsotope = IsotopeData("test", self.test_data)
        testIsotope.perform_deadtime_correction(dwell_time = self.dwell_time,
                                                dead_time = self.dead_time)
        testIsotope.perform_deadtime_correction(dwell_time = self.dwell_time,
                                                dead_time = self.dead_time)

    def test_mask_lt(self):
        testIsotope = IsotopeData("test", self.test_data)
        testIsotope.perform_deadtime_correction(dwell_time = self.dwell_time,
                                                dead_time = self.dead_time)
        y = testIsotope < 100;
        answer = np.array([[[False, False, False ],
                            [ True , False, False],
                            [ True,   True,  True]],

                           [[  True,  True,  True ],
                            [ True ,  True,  True],
                            [False, False , False]]])
        assert_true(np.array_equal(y, answer))

    def test_mask_gt(self):
        testIsotope = IsotopeData("test", self.test_data)
        testIsotope.perform_deadtime_correction(dwell_time = self.dwell_time,
                                                dead_time = self.dead_time)
        y = testIsotope > 100;
        answer = np.array([[[True, True, True ],
                            [ False , True, True],
                            [ False,   False,  False]],

                           [[  False,  False,  False ],
                            [ False ,  False,  False],
                            [True, True , True]]])
        assert_true(np.array_equal(y, answer))
        
        
    def test_importer_deadtime(self):
        testIsotope = IsotopeData("test", self.test_data)
        test_importer = Importer()
        test_importer.add_isotope(testIsotope)
        test_importer.deadtime_correct_all(dwell_time = self.dwell_time,
                                                  dead_time = self.dead_time)
                                                  
        assert_true(np.allclose(self.dt_corrected,
                                test_importer.get_isotope("test")._data))
