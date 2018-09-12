"""

.. module:: data_structures
    :synopsis: Data structure for Isotope and Ratio data.

.. moduleauthor:: Joshua Rehak <jsrehak@berkeley.edu>

"""

import numpy as np
import numpy.ma as ma
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class IsotopeData(object):
    """ Create an IsotopeData file for an isotope with given name, and data.
    
    :param isotope_label: Name of the isotope.
    :type isotope_label: string

    :param isotope_data: The data for the given isotope.
    :type isotope_data: 3D `numpy` array
    """

    def __init__(self, isotope_label, isotope_data):
        self._label = isotope_label
        self._data = np.array(isotope_data, dtype=float)
        self._is_deadtime_corrected = False

    def get_label(self):
        return self._label

    def __leq__(self, value):
        return self._data <= value
    
    def __lt__(self, value):
        return self._data < value

    def __gt__(self, value):
        return self._data > value

    def get_data(self):
        return self._data

    def get_mask(self, lower=0, upper=np.Inf):
        """Return a mask that will mask all data outside of the bounds given. \
           note that the numpy mask sets values that *will* be masked to True.

        :param lower: lower bound for mask (default 0), values less than this \
                      will be masked.
        :type type: float

        :param upper: upper bound for mask (default infinity) , values more \
                      than this will be masked.
        :type upper: float
        
        """
        return np.logical_or(self._data <= lower, self._data > upper)
    
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
        # Prevent applying deadtime correction twice
        if self._is_deadtime_corrected:
            raise RuntimeError("Error: Isotope " + self._label +
                               " is already deadtime corrected")

        print("Deadtime correction: isotope: " + self._label + ";\tdwell_time: " +
              str(dwell_time) + ";\tdead_time: " + str(dead_time))

        # Perform deadtime correction
        count_rate = np.divide(self._data, dwell_time)
        self._data = np.divide(count_rate,
                               1 - np.multiply(count_rate, dead_time))
        self._data *= dwell_time
        self._is_deadtime_corrected = True

    def plot(self, mask=None):
        """ Plot the isotope, with desired mask.

        :param mask: Mask to be used.
        :type mask: numpy bool array
        """

        # Get plot data
        if type(mask) is np.ndarray:
            plot_data = ma.masked_array(self._data, mask = mask)
        else:
            plot_data = self._data
        
        x = range(np.shape(plot_data)[1])
        y = range(np.shape(plot_data)[2])
        z_max = np.shape(plot_data)[0]
        X, Y = np.meshgrid(x,y)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')
        
        vmin = np.min(plot_data)
        vmax = np.max(plot_data)
        
        for z in range(z_max):
            x_y_data = plot_data[z, :, :]
            contour = ax.contourf(X, Y, x_y_data, zdir='z', offset=-z, alpha = 1,
                                  cmap = "CMRmap", vmin=vmin, vmax=vmax)
            
        cbar = fig.colorbar(contour, ax = ax)
        
        cbar.ax.set_xlabel("Counts")

        ax.set_zlim((-z_max,0))
        ax.set_zlim((-z_max,0))
        ax.set_xlabel("X (px)")
        ax.set_ylabel("Y (px)")
        ax.set_zlabel("Planes")
        ax.set_aspect('equal')
        ax.invert_xaxis()
        ax.view_init(elev=45, azim=45)

        values = "Min value: " + str(vmin) + " Max value: " + str(vmax);
        fig.text(x=0.2, y=0.02, s=values)
        
        plt.show()

    def sum(self, mask=None):
        """ Returns the sum of all the data in the dataset, with optional masking.

        :param mask: numpy mask array.
        :type mask: numpy array, optional)
        """
        masked_array = np.ma.array(self._data, mask=mask)
        return masked_array.sum()
        
    def __str__(self):
        return_string = "label: " + self._label + "; "
        return_string += "\tData size: " + str(np.shape(self._data))
        return_string += "\n\t Corrections: "
        if self._is_deadtime_corrected:
            return_string += "deadtime"
        return return_string

class RatioData(IsotopeData):
    def __init__(self, label, numerator_isotope, denomenator_isotope):
        self._label = label
        numerator_data = numerator_isotope.get_data()
        denomenator_data = denomenator_isotope.get_data()
        self._data = np.divide(numerator_data, denomenator_data,
                               out=np.zeros_like(denomenator_data),
                               where=denomenator_data!=0)
    
    def perform_deadtime_correction(self, dwell_time, dead_time):
        raise RuntimeError("Deadtime correction cannot be performed on a ratio")
