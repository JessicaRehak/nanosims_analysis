#!/usr/bin/env python3

#####################
# IDP ANALYSIS SCRIPT
#
# This shows data reduction for a non-standard (unknown) NanoSIMS analysis.
# File should be in the same directory where this script is.
#

# Import the required objects

import numpy as np
from nanosims_analysis.importer import Importer
from nanosims_analysis.data_structures import RatioData

# Set the filename
filename  = "./chim02 FIB SC Olivine_2.im"

print("Masked IDP analysis for file: " + str(filename))

# Create importer to import data, and import the data
importer = Importer()
importer.import_file(filename)

# Perform deadtime correction on all isotope data
importer.deadtime_correct_all(dead_time = 44*10**-9)

# Trim datasets if needed
trim = input("Trim front of dataset? [y/n] ")
if str(trim) == 'y':
    trim_amount = input("Please input number of cycles to trim: ")
    importer.trim_front_all(trim_amount)

# # Assign IsotopeData objects to each isotope
O16  = importer.get_isotope("16O")
O17  = importer.get_isotope("17O")
O18  = importer.get_isotope("18O")
Si28 = importer.get_isotope("28Si")
S32  = importer.get_isotope("32S")
Mg24 = importer.get_isotope("24Mg 16O")


    
# Makes a RatioData object
O18_to_O16 = RatioData("O18 to O16", numerator_isotope = O18,
                       denominator_isotope = O16)
O17_to_O16 = RatioData("O17 to O16", numerator_isotope = O17,
                       denominator_isotope = O16)


# Output to a VTK File.
#
# For this dataset, there is one pixel at the END of the x-range and y-range
# that need to be moved to the front of the data set (this is an issue with the
# .im file), so we set x_roll and y_roll to 1. Negative values will move that
# number of pixels from the front to the end, the reverse operation.
#
O16.to_VTK("./structured", x_roll = 1, y_roll = 1)
