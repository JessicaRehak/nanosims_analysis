"""

.. module:: data_structures
    :synopsis: Data structure for Isotope and Ratio data.

.. moduleauthor:: Joshua Rehak <jsrehak@berkeley.edu>

"""

class IsotopeData(object):
    """ Create an IsotopeData file for an isotope with given name, and data.
    
    :param isotope_name: Name of the isotope.
    :type isotope_name: string

    :param isotope_data: The data for the given isotope.
    :type isotope_data: 3D `numpy` array
    """

    def __init__(self, isotope_name, isotope_data):
        self._isotope_name = isotope_name
        self._isotope_data = isotope_data

    def __str__(self):
        return_string = "Isotope: " + self._isotope_name + "; "
        return_string += "Data size: " + str(np.shape(self._isotope_data))
        return return_string        
