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

    :param deadtime_correction: If True (default), performs deadtime correction on all
                     imported data.
    :type deadtime_correction: bool
    """
    def __init__(self, filename, deadtime_correction = True):
        self._isotope_data = {}
        self._filename = filename
        
        # Import data
        self._import_data()

    def _import_data(self):
        # Verify path exists
        if not Path(self._filename).is_file():
            raise RuntimeError('Bad filename')

        self._sims_object = sims.SIMS(self._filename)

        # Create an IsotopeData object for each isotope in the file
        for i, isotope_data in enumerate(self._sims_object.data):
            label = self._sims_object.header["label list"][i]
                      
            self._isotope_data.update({
                label:
                IsotopeData(label, isotope_data)})
