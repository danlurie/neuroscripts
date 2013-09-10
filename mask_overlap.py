"""
Python script to calculate the overlap between a group mask and 
individual subject masks. Useful for identifying outliers to be
excluded when creating group masks for whole-brain analyses.

I usually use a 90 percent group mask, and exclude any subjects with an overlap
of less than 95%.

Usage:
python mask_overlap.py /path/to/mask_list.txt /path/to/group_mask.nii.gz output_name

mask_list.txt is a file listing paths to individual subject masks, one per line

group_mask.nii.gz is usually a 90 percent group mask

output_name is the basename of the output file (e.g. "foo" will result in foo.csv)

Written by Dan Lurie (danjlurie@gmail.com)
"""

import numpy as np
import nibabel
import sys
import re

# Read in the command line arguments
a = sys.argv[1]
b = sys.argv[2]
output_name = sys.argv[3]

# Open a text file containing paths to individual masks
# and tell Python to be nice about line-endings
subject_list = open(a, "U")

# Load the group mask file
group_mask_file = nibabel.load(b)

# Open the group mask file
group_mask_image = group_mask_file.get_data()

# If group mask is not float 32, make it so
if group_mask_image.dtype != "float32": 
	group_mask_image = group_mask_image.astype(np.float32)

# Convert the group mask image to a 1D array
group_mask_1D = group_mask_image.flatten()

# Count the number of voxels within the group mask
group_mask_count = np.sum(group_mask_1D)

# Create a new file in which to write our output and tell
# Python to append new lines instead of overwriting the file
overlap_table = open("%s.csv" % output_name,"a")

for subject in subject_list:

	subject = subject.rstrip('\n')

	# Load the subject mask file
	subject_mask_file = nibabel.load(subject)

	# Open the group mask file
	subject_mask_image = subject_mask_file.get_data()

	# If subject mask is not float 32, make it so
	if subject_mask_image.dtype != "float32": 
		subject_mask_image = subject_mask_image.astype(np.float32)

	# Convert the group mask image to a 1D array
	subject_mask_1D = subject_mask_image.flatten()

	# Multiply the group mask by the subject mask
	overlap_array = np.dot(group_mask_1D,subject_mask_1D)

	# Count the number of overlapping voxels
	overlap_count = np.sum(overlap_array)

	# Calculate the percent overlap [do I need to convert to floats?]
	percent_overlap = overlap_count / group_mask_count

	# Write Subject ID and percent overlap to .csv file
	overlap_table.write("%s,%s\n" % (subject, percent_overlap))

# Close any files we opened, so they don't hog resources
overlap_table.close()
subject_list.close()






