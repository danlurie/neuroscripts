"""
Python script to count the number of subjects exhibiting each DSM Diagnosis.

To be run on data tables exported from the COINS system (http://ni.mrn.org)

Output is coins_diagnoses.csv

Usage:
python dx_count.py /path/to/dx_list.txt /path/to/coins_data.csv

Written by Dan Lurie (danjlurie@gmail.com)
"""

import rpy2.robjects as robjects
import sys
import re

# Read in the command line arguments
list_file = sys.argv[1]
data_file = sys.argv[2]

# Read in a text file containing diagnoses, and tell python to be nice about line endings
dx_list = open("%s" % list_file, "U")

# Read the COINS data csv into R
dx_table_raw = robjects.r("read.csv('%s')" % data_file)

# Create a python object of the R dataframe
dx_table = dx_table_raw.r_repr()

# Open a file where we will write script output
output_file = open("coins_diagnoses.csv", "a'")

# Write column headers to the output file.
output_file.write("Diagnosis Name, Total")

# For each diagnosis...
for dx in dx_list:
    # Strip \n from the end of each diagnosis
    dx = dx.rstrip('\n')
    # List the subjects with that diagnosis
    dx_rows = robjects.r("apply(%s, 1, function(x) any(x=='%s'))" % (dx_table, dx))
    # Count the subjects using R
    dx_sum = robjects.r('length(which(%s))' %dx_rows.r_repr())
    # Get R output
    dx_sum_string = dx_sum.r_repr()
    # Strip cruft from R output
    dx_print = re.sub("[^0-9]", "", dx_sum_string)
    # Write the diagnosis and number of subjects to the output file
    output_file.write("%s,%s\n" % (dx, dx_print))

# Close the output file
output_file.close()
