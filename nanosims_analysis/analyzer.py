"""

.. module:: analyzer
    :synopsis: Main program object, holds all data.

.. moduleauthor:: Joshua Rehak <jsrehak@berkeley.edu>


"""
import numpy as np
from pathlib import Path
from tkinter import Tk
from tkinter.filedialog import askopenfilename
import sims

class Analyzer():
    """ An object that uploads and stores data from a NanoSIMS file and allows
    the user to perform analysis.

    :param filename: If empty (default), will open a file menu to select a file \
                     to open. Otherwise, attempts to open this filename.
    :type filename: string, optional
    """
    def __init__(self, filename=""):
        # Initialize empty isotope data
        self._isotope_data = {}
        self._import_data(filename)

    def print_data(self):
        # Prints information about the data files stored in the Analyzer object
        for label, isotope in self._isotope_data.items():
            print(isotope)
        
    def _import_data(self, filename):
        # Get or save filename
        if not filename:
            root = Tk()
            root.withdraw()
            filename = askopenfilename()

        self._filename = filename
            
        # Verify file exists
        file_path = Path(self._filename)
        if not file_path.is_file():
            raise RuntimeError('Bad filename')

        # Create SIMS object, get all the data
        self._sims_object = sims.SIMS(self._filename)

        # Create an IsotopeData object for each isotope in the file
        for i, isotope_data in enumerate(self._sims_object.data):
            label = self._sims_object.header["label list"][i]
                      
            self._isotope_data.update({
                label:
                IsotopeData(label, isotope_data)})

class IsotopeData():
    """ An object containing all the data for a single isotope
    """

    def __init__(self, isotope_name, isotope_data):
        self._isotope_name = isotope_name
        self._isotope_data = isotope_data

    def __str__(self):
        return_string = "Isotope: " + self._isotope_name + "; "
        return_string += "Data size: " + str(np.shape(self._isotope_data))
        return return_string

    
