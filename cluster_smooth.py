"""
Python script to run spatial smoothing on multiple subjects using an SGE server.

Usage:
python cluster_smootb.py /path/to/subject_list.txt /path/to/input_dir sigma /path/to/mask.nii.gz /path/to/output_dir

subject_list.txt is a file listing subject IDs, one per line

input_dir is the directory containing unsmoothed functional files

sigma is the sigma of the guassian kernel to be used

mask.nii.gz is the group mask to be used

output_dir is the directory where smoothed files will be written

Written by Dan Lurie (danjlurie@gmail.com)
"""

import sys
import re
# Read in the command line arguments
a = sys.argv[1]
input_dir = sys.argv[2]
sigma = sys.argv[3]
mask_path = sys.argv[4]
output_dir = sys.argv[5]

# Open a text file containing subject IDs
# and tell Python to be nice about line-endings
subject_list = open(a, "U")

for sub in subject_list:
	sub = sub.rstrip('\n')
	commands = "fslmaths %s/%s.nii.gz -kernel -gauss %s -fmean -mas %s %s/%s_smoothed.nii.gz" % (input_dir, sub, sigma, mask_path, output_dir, sub)
	command_script = "qsub/scripts/%s.bash" % sub 
	command_script_file = open(command_script, 'w')
	command_script_file.write("#!/usr/bash\n")
	command_script_file.write(commands + "\n")
	command_script_file.close()
	log_file = "qsub/logs/%s.log" % sub
	queue_command = "qsub -S /bin/bash -V -cwd -o %s -j y %s" % (log_file, command_script_file)
	print queue_command
	os.system(queue_command)