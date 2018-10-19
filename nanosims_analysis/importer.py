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
    """ Importer object for importing data from a NanoSIMS file.
    """
    def __init__(self):
        self._isotopes = {}

    def add_isotope(self, isotope_data):
        """ Directly add an IsotopeData object to the importer

        :param isotope_data: IsotopeData object to add.
        :type isotope_data: IsotopeData
        """
        self._isotopes.update({
            isotope_data.get_label() : isotope_data})
    
    def get_isotope(self, label):
        """ Get the isotope identified by label

        :param label: isotope label
        :type label: string
        """
        return self._isotopes[label]

    def import_file(self, filename):
        """ Uploads and stores data from a NanoSIMS file.
        
        :param filename: Attempts to open this file to import data.
        :type filename: string

        """
        self._filename = filename
        
        # Verify path exists
        if not Path(self._filename).is_file():
            raise RuntimeError('Bad filename')

        self._sims_object = sims.SIMS(self._filename)

        # Create an IsotopeData object for each isotope in the file
        for i, isotope_data in enumerate(self._sims_object.data):
            label = self._sims_object.header["label list"][i]
                      
            self._isotopes.update({
                label:
                IsotopeData(isotope_label = label,
                            isotope_data = isotope_data)})

        
    def deadtime_correct_all(self, dead_time, dwell_time=0):
        """ Performs deadtime correction on each data set. Dwell time can usually be \
        found in the header for the NanoSIMS file but dead time must be given. \
        See :meth:`~nanosims_analysis.data_structures.IsotopeData.perform_deadtime_correction`.

        :param dead_time: Dead time in **seconds**.
        :type dead_time: float
        """
        if not dwell_time:
            self._dwell_time = float(
                self._sims_object.header["BFields"][0]["time per pixel"])
        else:
            self._dwell_time = dwell_time
        self._dead_time = dead_time

        for label, isotope in self._isotopes.items():
            isotope.perform_deadtime_correction(dwell_time = self._dwell_time,
                                                dead_time = self._dead_time)

    def trim_front_all(self, n):
        for label, isotope in self._isotopes.items():
            isotope.trim_front(int(n))

    def __str__(self):
        return_string = "Importer object\nImported file: " + self._filename + "\n";
        return_string += "Isotopes:\n"
        for label, data in self._isotopes.items():
            return_string += str(data) + "\n"
        return return_string[:-1]
