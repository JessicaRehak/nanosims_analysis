#!/usr/bin/env python3

#####################
# EXAMPLE SCRIPT
#
# This shows data reduction for a file called "chim02 FIB SC Olivine_2.im" in
# the directory where this file is.
#

# Import the required objects
from nanosims_analysis.importer import Importer
from nanosims_analysis.data_structures import RatioData

# Set the filename
filename  = "./chim02 FIB SC Olivine_2.im"

# Create importer to import data, and import the data
importer = Importer()
importer.import_file(filename)

# Perform deadtime correction on all isotope data
importer.deadtime_correct_all(dead_time = 44*10**-9)

# # Assign IsotopeData objects to each isotope
O16  = importer.get_isotope("16O")
O17  = importer.get_isotope("17O")
O18  = importer.get_isotope("18O")
Si28 = importer.get_isotope("28Si")
S32  = importer.get_isotope("32S")
Mg24 = importer.get_isotope("24Mg 16O")

# Makes a RatioData object
O18_to_O16 = RatioData("O18 to O16", numerator_isotope = O17,
                       denominator_isotope = O16)

# Calculate R17 = O17/O16 using sums
# Generate mask from O16 data
O16mask = O16.get_mask(lower = 500)
# Calculate the ratio using the masked sums
R17 = O17.sum(O16mask)/O16.sum(O16mask)
