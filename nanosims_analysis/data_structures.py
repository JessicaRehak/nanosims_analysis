"""

.. module:: data_structures
    :synopsis: Data structure for Isotope and Ratio data.

.. moduleauthor:: Joshua Rehak <jsrehak@berkeley.edu>

"""

import numpy as np

class IsotopeData(object):
    """ Create an IsotopeData file for an isotope with given name, and data.
    
    :param isotope_label: Name of the isotope.
    :type isotope_label: string

    :param isotope_data: The data for the given isotope.
    :type isotope_data: 3D `numpy` array
    """

    def __init__(self, isotope_label, isotope_data):
        self._label = isotope_label
        self._data = isotope_data

    def perform_deadtime_correction(self, dwell_time, dead_time):
        r""" Perform deadtime correction on the count data:\

        .. math:: n=\frac{n_0}{1-n_0\tau}
        
        using dead time, :math:`\tau`, and :math:`n_0` is in counts per second:\

        .. math:: n_0=\frac{n}{T}

        where :math:`n` is the counts, and :math:`T` is the dwell time.

        .. important::  Dwell time and dead time may be in any unit of time as \
           long as they are in the *same* units.
        
        :param dwell_time: Dwell time.
        :type dwell_time: float

        :param dead_time: Dead time.
        :type dead_time: float
        """
        # Get count rates
        count_rate = np.divide(self._data, dwell_time)
        self._data = np.divide(count_rate,
                               1 - np.multiply(count_rate, dead_time))
        self._label += " (deadtime corrected)"
        
    def __str__(self):
        return_string = "label: " + self._label + "; "
        return_string += "Data size: " + str(np.shape(self._data))
        return return_string
