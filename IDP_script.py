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
filename  = "./Chim03 HumptyDumpty-1.im"

print("Masked IDP analysis for file: " + str(filename))

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
O18_to_O16 = RatioData("O18 to O16", numerator_isotope = O18,
                       denominator_isotope = O16)
O17_to_O16 = RatioData("O17 to O16", numerator_isotope = O17,
                       denominator_isotope = O16)


# TO DO: Correct for QSA on RatioData object 
# TO DO: Add reduced resolution option for imaging

# Generate mask from O16 data
O16mask = O16.get_mask(lower = 1000)

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
primary_current = int(input("Please input primary current (in pA): "))
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

d17OstdM = 51.37322054
d17OstdErr = 2.729224551
d18OstdM = 78.87885747
d18OstdErr = 1.207378055

d17OstdL = 2.7             # San Carlos olivine (Tanaka and Nakamura, 2013) 
d18OstdL = 5.3             # San Carlos olivine (Tanaka and Nakamura, 2013) 
print("Standardized to San Carlos olivine (Tanaka and Nakamura, 2013)")

# TO DO: Add python prompt for standard inputs above instead of being hard-written into script
# TO DO: Add uncertainties for measured standards (use standard deviation)

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
