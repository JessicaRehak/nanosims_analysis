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
#from src.hl import gridToVTK 

# Set the filename
filename  = "Chim05_SC_Olivine_FIB_2_1.im"

print("Masked IDP analysis for file: " + str(filename))

# Create importer to import data, and import the data
importer = Importer()
importer.import_file(filename)

# Ask user for primary current
primary_current = int(input("Please input primary current (in pA): "))

# Perform deadtime correction on all isotope data
importer.deadtime_correct_all(dead_time = 44*10**-9)

# For this dataset, there is one pixel at the END of the x-range and y-range
# that need to be moved to the front of the data set (this is an issue with the
# .im file), so we set x_roll and y_roll to 1. Negative values will move that
# number of pixels from the front to the end, the reverse operation.
importer.roll_all(x_roll = 1, y_roll = 1)

# # Assign IsotopeData objects to each isotope
O16  = importer.get_isotope("16O")
O17  = importer.get_isotope("17O")
O18  = importer.get_isotope("18O")
Si28 = importer.get_isotope("28Si")
S32  = importer.get_isotope("32S")
Mg24 = importer.get_isotope("24Mg 16O")

# Trim datasets
trim = input("Trim front of dataset? [y/n] ")
if str(trim) == 'y':
    trim_amount = input("Please input number of cycles to trim from front: ")
    importer.trim_front_all(trim_amount)

trimb = input("Trim back of dataset? [y/n] ")
if str(trimb) == 'y':
    trimb_amount = input("Please input number of cycles to trim from back: ")
    importer.trim_back_all(trimb_amount)
    
# return number of cycles
ncycles = O16.n_cycles()
    
# Makes a RatioData object
O18_to_O16 = RatioData("O18 to O16", numerator_isotope = O18,
                       denominator_isotope = O16)
O17_to_O16 = RatioData("O17 to O16", numerator_isotope = O17,
                       denominator_isotope = O16)


# TO DO: Correct for QSA on RatioData object 
# TO DO: Add reduced resolution option for imaging

# Generate mask from O16 data
maskQ = input("Mask data in 16O? [y/n] ")
if str(maskQ) == 'y':
    mask_low = int(input("Please input lower limit on 16O (counts/pixel): "))
else:
    mask_low = 0
    
O16mask = O16.get_mask(lower = mask_low)

# Return number of pixels in mask
pixels = O16.n_pixels(O16mask)

# Calculate total counts by sums
O16tot = O16.sum(O16mask)
O17tot = O17.sum(O16mask)
O18tot = O18.sum(O16mask)

# Calculate the ratio using the masked sums
R17init = O17tot/O16tot
R18init = O18tot/O16tot


# Correct for QSA on bulk Ratio using methods of Hillion et al., 2008
# Corrected Ratio = ratio_measured/(1+beta*K)
# where K = O16tot / (primary current * 6.2415e6 * dwell time * N pixels)   

dwell_time = 3000  # TO DO: edit so that it reads dwell time, instead of being hard wired into script
beta = 0.75
K = O16tot / (primary_current * 6.2415*10**6 * dwell_time * pixels)
R17 = R17init / (1 + beta*K)
R18 = R18init / (1 + beta*K)
print("QSA correction using beta value: " + str(beta) + " (Hillion et al., 2008)")

# Calculate uncertainties 
# Instrumental statistical uncertainty for isotope = root(N) where N is total counts
rootN16 = np.sqrt(O16tot)
rootN17 = np.sqrt(O17tot)
rootN18 = np.sqrt(O18tot)

# Statistical uncertainty for ratio
sigR17 = R17 * np.sqrt((rootN17/O17tot)**2 + (rootN16/O16tot)**2)
sigR18 = R18 * np.sqrt((rootN18/O18tot)**2 + (rootN16/O16tot)**2)

# For a more conservative estimate, use standard error -- stdev / sqrt(n) for Ratio data where n is number of analyses
stdev17 = np.std(O17_to_O16.get_data(O16mask))
stdev18 = np.std(O18_to_O16.get_data(O16mask))

sig17 = stdev17 / np.sqrt(pixels)
sig18 = stdev18 / np.sqrt(pixels)

# Quadratically combined Ratio error - maybe overcounts the error? Will probably not use this value
Qsig17 = np.sqrt(sig17**2 + sigR17**2)
Qsig18 = np.sqrt(sig18**2 + sigR18**2)

# define standard ratio values for VSMOW (Baertschi, 1976; Fahey et al., 1987)
R17smow = 0.00038288
R18smow = 0.0020052

# Calculate the delta values for ratio and uncertainties
delta17O = ((R17/R17smow)-1)*1000  # For bulk ratio
delta18O = ((R18/R18smow)-1)*1000  # For bulk ratio

RDatad17O = ((O17_to_O16.get_data(O16mask)/R17smow)-1)*1000  # For RatioData object
RDatad18O = ((O18_to_O16.get_data(O16mask)/R18smow)-1)*1000  # For RatioData object

dstdev17 = np.std(RDatad17O)
dstdev18 = np.std(RDatad18O)

dsig17 = (sig17/R17smow)*1000
dsig18 = (sig18/R18smow)*1000


d2sig17 = 2*dsig17
d2sig18 = 2*dsig18

# Correct for instrumental mass fractionation (IMF)
# Input standard (std) delta values, as measured by SIMS (M), and literature values (L) in units of permil

##d17OstdM = float(input("Please input measured delta 17O value for standard (in permil): "))
##d17OstdErr = float(input("Please input 1 sigma uncertainty on d17O for standard (in permil): "))
##d18OstdM = float(input("Please input measured delta 17O value for standard (in permil): "))
##d18OstdErr = float(input("Please input 1 sigma uncertainty on d18O for standard (in permil): "))

d17OstdM = 18.62939811
d17OstdErr = 3.102282541
d18OstdM = 42.72560561
d18OstdErr = 1.366749018

d17OstdL = 2.7             # San Carlos olivine (Tanaka and Nakamura, 2013) 
d18OstdL = 5.3             # San Carlos olivine (Tanaka and Nakamura, 2013) 
print("Standardized to San Carlos olivine (Tanaka and Nakamura, 2013)")

IMF17 = d17OstdM - d17OstdL
IMF18 = d18OstdM - d18OstdL

d17Ocorr = delta17O - IMF17
d18Ocorr = delta18O - IMF18

# Calculate combined uncertainties for IMF corrected data 
Err17 = np.sqrt(np.power(dsig17,2) + np.power(d17OstdErr,2))
Err18 = np.sqrt(np.power(dsig18,2) + np.power(d18OstdErr,2))

twoErr17 = 2 * Err17
twoErr18 = 2 * Err18

# Print results (Needs to be updated with correct uncertainties)
print("The IMF corrected d17O value is: " + str(d17Ocorr) + " +/- " + str(twoErr17) + " permil (2 sigma)" + "\nThe IMF corrected d18O value is: " + str(d18Ocorr) + " +/- " + str(twoErr18) + " permil (2 sigma)")

# TO DO: Append Data to CSV - write function that appends
# append = input("Append data to CSV file? [y/n] ")

# Save numpy array to file 
# from tempfile import TemporaryFile
# outfile = TemporaryFile()
# np.save("O16.npy", O16, fix_imports=True)


## Export data to VTK file 
# Dimensions 
O16.to_VTK("./VTK_output_files/" + filename + "_16O", x_roll = 1, y_roll = 1)

O18_to_O16.to_VTK("./VTK_output_files/" + filename + "_18to16O", x_roll = 1, y_roll = 1, mask=O16mask)

O16.plot(O16mask)




