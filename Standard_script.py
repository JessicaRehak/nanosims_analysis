#!/usr/bin/env python3

#####################
# STANDARD ANALYSIS SCRIPT
#
# This shows data reduction for standard analyses. File should be located in
# the directory where this file is. 
#

# Import the required objects

import numpy as np
from nanosims_analysis.importer import Importer
from nanosims_analysis.data_structures import RatioData

# Set the filename
filename  = "./im21 SC olivine FIB .im"

print("Standard analysis for file: " + str(filename))

# Ask for primary current
primary_current = int(input("Please input primary current (in pA): "))

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

# Trim datasets
trim = input("Trim front of dataset? [y/n] ")
if str(trim) == 'y':
    trim_amount = input("Please input number of cycles to trim: ")
    for dataset in [O16, O17, O18, Si28, S32, Mg24]:
        dataset.trim_front(int(trim_amount))

# Makes a RatioData object
O18_to_O16 = RatioData("O18 to O16", numerator_isotope = O18,
                       denominator_isotope = O16)
O17_to_O16 = RatioData("O17 to O16", numerator_isotope = O17,
                       denominator_isotope = O16)


# TO DO: Correct for QSA on RatioData object (point by point ratios)
#    if ratio_label in ['17O/16O', '18O/16O', '32S/33S', '32S/34S', '32S/36S' ]:
#        beta = 0.75
#    elif ratio_label in ['12C/13C', '12C2/13C12C']:
#        beta = 1
#    elif ratio_label in ['28Si/29Si', '28Si/30Si']:
#        beta =0.6
#    else:
#        # Theoretical value for Poisson stats
#        beta = 0.5
# Need to make QSA corrected RatioData object, where value is adjusted to
# beta = 0.75
# ratio_corr = ratio_measured/(1+beta*K)

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

RDatad17O = ((O17_to_O16.get_data(O16mask)/R17smow)-1)*1000  # delta values for RatioData object
RDatad18O = ((O18_to_O16.get_data(O16mask)/R18smow)-1)*1000  # delta values for RatioData object

dstdev17 = np.std(RDatad17O)  # standard deviation in delta values
dstdev18 = np.std(RDatad18O)  # standard deviation in delta values

dsig17 = (sig17/R17smow)*1000  # convert error to delta values
dsig18 = (sig18/R18smow)*1000  # convert error to delta values

d2sig17 = 2*dsig17  # 2 sigma
d2sig18 = 2*dsig18  # 2 sigma

# Print results 
print("The d17O value is: " + str(delta17O) + " +/- " + str(d2sig17) + " permil (2 sigma)" + "\nThe d18O value is: " + str(delta18O) + " +/- " + str(d2sig18) + " permil (2 sigma)")

O16.plot(O16mask)