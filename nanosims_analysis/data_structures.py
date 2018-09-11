"""

.. module:: data_structures
    :synopsis: Data structure for Isotope and Ratio data.

.. moduleauthor:: Joshua Rehak <jsrehak@berkeley.edu>

"""

import numpy as np
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
        self._data = isotope_data
        self._is_deadtime_corrected = False

    def get_label(self):
        return self._label

    def __lt__(self, value):
        return self._data < value

    def __gt__(self, value):
        return self._data > value
    
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
        self._is_deadtime_corrected = True

    def plot(self):
        x = range(np.shape(self._data)[1])
        y = range(np.shape(self._data)[2])
        z_max = np.shape(self._data)[0]
        X, Y = np.meshgrid(x,y)
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')
        
        vmin = np.min(self._data)
        vmax = np.max(self._data)
        
        for z in range(z_max):
            x_y_data = self._data[z, :, :]
            contour = ax.contourf(X, Y, x_y_data, zdir='z', offset=-z, alpha = 1,
                                  cmap = "CMRmap", vmin=vmin, vmax=vmax)
            
        cbar = figp.colorbar(contour, ax = ax)
        
        if self._is_deadtime_corrected:
            cbar.ax.set_xlabel("Counts/sec")
        else:
            cbar.ax.set_xlabel("Counts")

        ax.set_zlim((-z_max,0))
        ax.set_zlim((-z_max,0))
        ax.set_xlabel("X (px)")
        ax.set_ylabel("Y (px)")
        ax.set_zlabel("Planes")
        ax.set_aspect('equal')
        ax.invert_xaxis()
        ax.view_init(elev=45, azim=45)
        
        plt.show()
        
    def __str__(self):
        return_string = "label: " + self._label + "; "
        return_string += "\tData size: " + str(np.shape(self._data))
        return_string += "\n\t Corrections: "
        if self._is_deadtime_corrected:
            return_string += "deadtime"
        return return_string
