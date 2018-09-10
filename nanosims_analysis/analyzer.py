"""

.. module:: analyzer
    :synopsis: Main program object, holds all data.

.. moduleauthor:: Joshua Rehak <jsrehak@berkeley.edu>

"""
from pathlib import Path
from tkinter import Tk
from tkinter.filedialog import askopenfilename
import isotope_data
import sims

class Analyzer():
    """ An object that uploads and stores data from a NanoSIMS file and allows
    the user to perform analysis.
    """
    def __init__(self):
        self.isotope_data = []
        self.filename = ""
        menu()
        
    def menu(self):
        print("NanoSIMS Analyzer")


        
        print("1. Open file")

    
    def __open_file(self):
        # Get filename
        if filename:
            self.filename = filename
        else:
            root = Tk()
            root.withdraw()
            self.filename = askopenfilename()
            
        # Verify file exists
        file_path = Path(self.filename)
        if not file_path.is_file():
            raise RuntimeError('Bad filename')

        # Create SIMS object, get all the data
        self.sims_object = sims.SIMS(self.filename)

        # Create an IsotopeData object for each isotope in the file
        for i, isotope_data in enumerate(self.sims_object.data):
            self.isotope_data.append(isotope_data.IsotopeData(
                self.sims_object.header["label list"][i],
                isotope_data))    

class IsotopeData():
    """ An object containing all the data for a single isotope
    """

    def __init__(self, isotope_name, isotope_data):
        self.isotope_name = isotope_name
        self.isotope_data = isotope_data
