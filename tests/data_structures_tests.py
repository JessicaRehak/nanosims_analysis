from nose.tools import *
import numpy as np
import numpy.ma as ma

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

        cls.dt_corrected = np.array([[[0.92673054, 0.98710203, 0.69649064],
                                      [0.13295208, 0.66058424, 0.71824281],
                                      [0.39265746, 0.03355262, 0.43079691]],

                                     [[0.0554341 , 0.59302009, 0.4309345 ],
                                      [0.46924814, 0.50525276, 0.38235408],
                                      [0.85084899, 0.63285645, 0.93021143]]])


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

    @raises(RuntimeError)
    def test_isotope_data_trim_too_much(self):
        for n in [4, 5, 100, 32]:
            random_data = np.random.rand(n, 20, 20)
            testIsotope = IsotopeData("test", random_data)
            testIsotope.trim_front(n+1)

    def test_isotope_data_trim_ok(self):
        for n in [4, 5, 100, 32]:
            random_data = np.random.rand(n, 20, 20)
            testIsotope = IsotopeData("test", random_data)
            t = np.random.randint(0, n-1)
            testIsotope.trim_front(t)
            assert_equal(np.shape(testIsotope._data)[0],
                         n - t)

    def test_get_mask(self):
        testIsotope = IsotopeData("test", self.test_data)
        testIsotope.perform_deadtime_correction(dwell_time = self.dwell_time,
                                                dead_time = self.dead_time)
        y = testIsotope.get_mask(lower = 0.1, upper = 0.7)
        answer = np.array([[[False, False, True ],
                            [ True , True, False],
                            [ True,   False,  True]],
                           
                           [[  False,  True,  True ],
                            [ True ,  True,  True],
                           [False, True , False]]])
        assert_true(np.array_equal(y, np.invert(answer)))

    def test_n_pixels(self):
        testIsotope = IsotopeData("test", self.test_data)
        assert_equal(testIsotope.n_pixels(), 18)

    def test_n_pixels_with_mask(self):
        testIsotope = IsotopeData("test", self.test_data)
        y = testIsotope.get_mask(lower = 0.1, upper = 0.7)
        assert_equal(testIsotope.n_pixels(mask = y), 11)
        
    def test_sum(self):
        testIsotope = IsotopeData("test", self.test_data)
        testIsotope.perform_deadtime_correction(dwell_time = self.dwell_time,
                                                dead_time = self.dead_time)
        sum = 9.82926987
        assert_true(np.isclose(sum, testIsotope.sum()))

    def test_sum_mask(self):
        testIsotope = IsotopeData("test", self.test_data)
        testIsotope.perform_deadtime_correction(dwell_time = self.dwell_time,
                                                dead_time = self.dead_time)
        y = testIsotope.get_mask(lower = 0.1, upper = 0.7)
        sum = 5.32714735
        assert_true(np.isclose(sum, testIsotope.sum(mask = y)))
        
    def test_lt(self):
        testIsotope = IsotopeData("test", self.test_data)
        testIsotope.perform_deadtime_correction(dwell_time = self.dwell_time,
                                                dead_time = self.dead_time)
        y = testIsotope < 0.5;
        answer = np.array([[[False, False, False ],
                            [ True , False, False],
                            [ True,   True,  True]],

                           [[  True,  False,  True ],
                            [ True ,  False,  True],
                            [False, False , False]]])
        assert_true(np.array_equal(y, answer))

    def test_gt(self):
        testIsotope = IsotopeData("test", self.test_data)
        testIsotope.perform_deadtime_correction(dwell_time = self.dwell_time,
                                                dead_time = self.dead_time)
        y = testIsotope > 0.5;
        answer = np.array([[[True, True, True ],
                            [ False , True, True],
                            [ False,   False,  False]],

                           [[  False,  True,  False ],
                            [ False ,  True,  False],
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

    def test_get_data(self):
        testIsotope = IsotopeData("test", self.test_data)
        testIsotope.perform_deadtime_correction(dwell_time = self.dwell_time,
                                                dead_time = self.dead_time)
        y = testIsotope.get_mask(lower = 0.1, upper = 0.7)
        answer = np.array([[[False, False, True ],
                            [ True , True, False],
                            [ True,   False,  True]],
                           
                           [[  False,  True,  True ],
                            [ True ,  True,  True],
                           [False, True , False]]])
        assert_true(np.allclose(testIsotope.get_data(), self.dt_corrected))
        assert_true(ma.allclose(testIsotope.get_data(mask = y),
                                ma.array(self.dt_corrected, mask = y)))
