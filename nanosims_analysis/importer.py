"""

.. module:: importer
    :synopsis: Imports data from a NanoSims file.

.. moduleauthor:: Joshua Rehak <jsrehak@berkeley.edu>

"""
from nanosims_analysis.data_structures import IsotopeData
import numpy as np
from pathlib import Path
import sims

class Importer(object):
    """ Uploads and stores data from a NanoSIMS file and performs deadtime \
    correction.

    :param filename: Attempts to open this file to import data.
    :type filename: string

    """
    def __init__(self, filename):
        self._isotopes = {}
        self._filename = filename
        
        # Import data
        self._import_data()

    def get_isotope(self, label):
        """ Get the isotope identified by label

        :param label: isotope label
        :type label: string
        """
        return self._isotopes[label]
        
    def perform_deadtime_correction(self, dead_time):
        """ Performs deadtime correction on each data set. Dwell time can be \
        found in the header for the NanoSIMS file but dead time must be given.

        :param dead_time: Dead time in **seconds**.
        :type dead_time: float
        """
        self._dwell_time = float(
            self._sims_object.header["BFields"][0]["time per pixel"])
        self._dead_time = dead_time

        for label, isotope in self._isotopes.items():
            isotope.perform_deadtime_correction(self._dwell_time, self._dead_time)

    def __str__(self):
        return_string = "Importer object\nImported file: " + self._filename + "\n";
        return_string += "Isotopes:\n"
        for label, data in self._isotopes.items():
            return_string += "  " + str(data) + "\n"
        return return_string[:-1]

    def _import_data(self):
        # Verify path exists
        if not Path(self._filename).is_file():
            raise RuntimeError('Bad filename')

        self._sims_object = sims.SIMS(self._filename)

        # Create an IsotopeData object for each isotope in the file
        for i, isotope_data in enumerate(self._sims_object.data):
            label = self._sims_object.header["label list"][i]
                      
            self._isotopes.update({
                label:
                IsotopeData(label, isotope_data)})
