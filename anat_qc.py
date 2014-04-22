'''
Program to calculate a number of quality control measures for anatomical neuroimaging data (T1 images).

Set up to run on CPAC (http://fcp-indi.github.io) output directories, but should be easily modifiable to
run on any data structure.

Output is written to a CSV, currently specified as a path on line 357 (this will be changed to a
command-line argument in future versions)

This script can be run either serially (on a single core, one subject at a time) or in parallel
(on multiple cores, multiple subjects at a time). By default, anat_qc.py will run on 10 cores.
To change the number of cores or switch to single-core mode, comment out the relevant code near
the bottom of the script.

Usage:

python anat_qc.py pipeline_path csf_thresh wm_thresh gm_thresh

pipeline_path is the path to the CPAC pipeline on which you would like to run the QC.
By default, all subejcts will be run (but you can also manually specifiy a subject list, see lines 282 and 283)
Note: When specifying pipeline_path, make sure it has a trailing slash

csf, wm, and gm_thresh are the tissue probability thresholds (e.g. CPAC strategy) to run on.
These should be specified as floats, just like in the CPAC pipeline configuration file.

Created by Dan Lurie (danjlurie@gmail.com) based on code from Cameron Craddock (algorithms/calculations)
'''


import numpy as np
import nibabel as nb
import pandas as pd
import scipy.ndimage as nd
import scipy.stats as stats
import os
import sys
from multiprocessing import Pool
import code
from functools import partial

class do_qc():

	def __init__(self, sub_id, pipeline_path, csf_threshold, gm_threshold, wm_threshold):
		self.sub_id = sub_id
		self.pipeline_path = pipeline_path
		self.csf_threshold = csf_threshold
		self.gm_threshold = gm_threshold
		self.wm_threshold = wm_threshold

		self.anat_path = self.pipeline_path + self.sub_id + '/anatomical_reorient/anat_resample.nii.gz'
		self.fg_mask_path = self.pipeline_path + self.sub_id + '/anatomical_reorient/qc_head_mask.nii.gz'
		self.csf_mask_path = self.pipeline_path + self.sub_id + '/anatomical_csf_mask/_csf_threshold_' + self.csf_threshold + '/segment_prob_0_maths_maths_maths.nii.gz'
		self.gm_mask_path = self.pipeline_path + self.sub_id + '/anatomical_gm_mask/_gm_threshold_' + self.gm_threshold + '/segment_prob_1_maths_maths_maths.nii.gz'
		self.wm_mask_path = self.pipeline_path + self.sub_id + '/anatomical_wm_mask/_wm_threshold_' + self.wm_threshold + '/segment_prob_2_maths_maths_maths.nii.gz'

		# Read in and load anatomical data.
		self.anat_image = self.load_image(self.anat_path, 'anatomical')
		self.anat_data = self.anat_image.get_data()

		# Create a head mask based on the anatomical data (but check first if this has already been created).
		if os.path.isfile(self.pipeline_path + self.sub_id + '/anatomical_reorient/qc_head_mask.nii.gz') == False:
			mask_command = '3dAutomask -q -prefix ' + self.pipeline_path + self.sub_id + '/anatomical_reorient/qc_head_mask.nii.gz -dilate 2 ' + self.pipeline_path + self.sub_id + '/anatomical_reorient/anat_resample.nii.gz'
			os.system(mask_command)
		else:
			pass

		# Read in, verify, and load foreground mask data.
		self.fg_mask_image = self.load_image(self.fg_mask_path, 'foreground mask')
		self.verify_mask(self.fg_mask_image, self.anat_image)
		self.fg_mask_data = self.fg_mask_image.get_data()

		# Read in, verify, and load GM mask data.
		self.gm_mask_image = self.load_image(self.gm_mask_path, 'GM mask')
		self.verify_mask(self.gm_mask_image, self.anat_image)
		self.gm_mask_data = self.gm_mask_image.get_data()

		# Read in, verify, and load WM mask data.
		self.wm_mask_image = self.load_image(self.wm_mask_path, 'WM mask')
		self.verify_mask(self.wm_mask_image, self.anat_image)
		self.wm_mask_data = self.wm_mask_image.get_data()

		# Read in, verify, and load CSF mask data.
		self.csf_mask_image = self.load_image(self.csf_mask_path, 'CSF mask')
		self.verify_mask(self.csf_mask_image, self.anat_image)
		self.csf_mask_data = self.csf_mask_image.get_data()

		# Initialize a dictionary to hold the measures we calcualte, and add sub_id as the first value
		self.qc_dict = dict()
		self.qc_dict['subject'] = self.sub_id

		# Calculate mean, standard deviation, and number of voxels for foreground.
		self.qc_dict['fg_mean'] = np.mean(self.anat_data[self.fg_mask_data == 1])
		self.qc_dict['fg_std'] = np.std(self.anat_data[self.fg_mask_data == 1])
		self.qc_dict['fg_size'] = self.fg_mask_data.sum()

		# Calculate mean, standard deviation, and number of voxels for background.
		self.qc_dict['bg_mean'] = np.mean(self.anat_data[self.fg_mask_data == 0])
		self.qc_dict['bg_std'] = np.mean(self.anat_data[self.fg_mask_data == 0])
		self.qc_dict['bg_size'] = self.fg_mask_data.size - self.fg_mask_data.sum()

		# Calculate mean, standard deviation, and number of voxels for WM.
		self.qc_dict['wm_mean'] = np.mean(self.anat_data[self.wm_mask_data == 1])
		self.qc_dict['wm_std'] = np.mean(self.anat_data[self.wm_mask_data == 1])
		self.qc_dict['wm_size'] = self.wm_mask_data.sum()

		# Calculate mean, standard deviation, and number of voxels for GM.
		self.qc_dict['gm_mean'] = np.mean(self.anat_data[self.gm_mask_data == 1])
		self.qc_dict['gm_std'] = np.mean(self.anat_data[self.gm_mask_data == 1])
		self.qc_dict['gm_size'] = self.gm_mask_data.sum()

		# Calculate mean, standard deviation, and number of voxels for CSF.
		self.qc_dict['csf_mean'] = np.mean(self.anat_data[self.csf_mask_data == 1])
		self.qc_dict['csf_std'] = np.mean(self.anat_data[self.csf_mask_data == 1])
		self.qc_dict['csf_size'] = self.csf_mask_data.sum()

		# Calculate SNR
		# SNR = (mean GM intensity) / (std of background intensities)
		self.qc_dict['snr'] = self.qc_dict['gm_mean'] / self.qc_dict['bg_std']

		# Calculate contrast to noise ratio (CNR)
		# CNR = (absolute difference between WM signal and GM signal) / (std of background intensities)
		self.qc_dict['cnr'] = np.abs(self.qc_dict['gm_mean'] - self.qc_dict['wm_mean']) / self.qc_dict['bg_std']

		# Calculate FBER
		self.qc_dict['fber'] = self.calc_fber(self.anat_data, self.fg_mask_data)

		# Calculate EFC
		self.qc_dict['efc'] = self.calc_efc(self.anat_data)

		# Calculate QI1 and handle exceptions if data is float.
		try: 
			self.qc_dict['qi1'] = self.calc_artifacts(self.anat_data, self.fg_mask_data)
		except:
			pass
		

	def load_image(self, image_file, image_type):
		# Try to read in an image file, give a helpful error message based on arguments if the read fails.
		try:
			image = nb.load(image_file)
		except IOError as e:
			print("Error reading %s file for subject %s:\n %s" % (image_type, self.sub_id, e.strerror))
			raise
		except:
			print("Error reading %s file for subject %s: UNKNOWN ERROR\n %s" % (image_type, self.sub_id, sys.exc_info()[0]))
			raise

		return image 

	def verify_mask(self, mask_image, anat_image):
		# Check that the specified mask is binary.
		mask_values = np.unique(mask_image.get_data())
		if (mask_values.size != 2) or not (mask_values == [0, 1]).all():
			print("Error: Mask for subject %s is not binary." % (self.sub_id))
			raise

		# Verify that the mask and anatomical images have the same dimensions.
		if anat_image.shape != mask_image.shape:
			print("Error: Mask and anatomical image for subject %s are different dimensions." % (self.sub_id))
			raise

		# Verify that the mask and anatomical images are in the same space (have the samme affine matrix)
		if (mask_image.get_affine() == anat_image.get_affine()).all == False:
			print("Error: Mask and anatomical image for subject %s are not in the same space" % (self.sub_id))
			raise

	def calc_fber(self, anat_data, mask_data):
		# Calculate Foreground:Background Energy Ratio
		
		# FBER = (mean foreground energy) / (mean background energy)
		fber = ((np.abs(anat_data[mask_data == 1]) ** 2).sum() / (mask_data.sum())) / \
		((np.abs(anat_data[mask_data == 0]) ** 2).sum() / (mask_data.size - mask_data.sum()))
		
		return fber

	def calc_efc(self, anat_data):
		# Calculate the Entropy Focus Criterion (Atkinson 1997, IEEE TMI)
		
		# We noramlize the original equation by the maximum entropy so our EFC can be 
		# easily compared across images with different dimensions.

		# Calculate the maximum value of the EFC (which occurs any time all voxels have the same value)
		efc_max = 1.0 * np.prod(anat_data.shape) * (1.0 / np.sqrt(np.prod(anat_data.shape))) * \
		np.log(1.0 / np.sqrt(np.prod(anat_data.shape)))

		# Calculate the total image energy
		b_max = np.sqrt((anat_data**2).sum())

		# Calculate EFC (add 1e-16 to the image data to keep log happy)
		efc = (1.0 / efc_max) * np.sum((anat_data / b_max) * np.log((anat_data + 1e-16) / b_max))

		return efc


	def calc_artifacts(self, anat_data, fg_mask_data):
		# Detect artifacts in the anatomical image using the method described in Mortamet et al. 2009 (MRM)
		# Calculates QI1, the fraction of total voxels that within artifacts.
		
		# Optionally, also calculates QI2, the distance between the distribution of noise voxel
		# (non-artifact background voxels) intensities, and a Ricean distribution.

		# Define the image background by taking the inverse of the foreground mask.
		bg_mask = (fg_mask_data == 0) * 1

		# Create an image containing only background voxels (everything outside bg_mask set to 0)
		background = anat_data.copy()
		background[bg_mask != 1] = 0


		# For now, our calculation can't handle floats, so we check that all values are integers.
		# (This is messy, need to find a better way to do this check.)
		if background.dtype == np.int or background.dtype == np.int16 or background.dtype == np.int32 \
		or background.dtype == np.int64:
			data_type = True

		else:
			data_type = False

		# If our data is of the right type, calculate QI1, otherwise throw an error.
		if data_type == True:

			# Find the background threshold (the most frequently occurring value excluding 0)
			bg_counts = np.bincount(background.flatten())
			bg_threshold = np.argmax(bg_counts[1:]) + 1

			# Apply this threshold to the background voxels to identify voxels contributing artifacts. 
			background[background <= bg_threshold] = 0
			background[background != 0] = 1

			# Create a structural element to be used in an opening operation.
			struct_elmnt = np.zeros((3,3,3))
			struct_elmnt[0,1,1] = 1
			struct_elmnt[1,1,:] = 1
			struct_elmnt[1,:,1] = 1
			struct_elmnt[2,1,1] = 1

			# Perform an opening operation on the background data.
			background = nd.binary_opening(background, structure=struct_elmnt)

			# Count the number of voxels that remain after the opening operation. These are artifacts.
			QI1 = background.sum() / float(bg_mask.sum())

			"""
			# Now lets focus on the noise, which is everything in the background that was not identified as artifact
			bgNoise = anatData[(bgMask-bg)==1]

			# calculate the histogram of the noise and its derivative
			H = np.bincount(bgNoise)
			H=1.0*H/H.sum()
			dH = H[1:]-H[:-1]
	
			# find the first value on the right tail, i.e. tail with negative slope, i.e. dH < 0 that is less than or equal to half of the histograms max
			firstNegSlope=np.nonzero(dH<0)[0][0]
			halfMaxRightTail=np.nonzero(H[firstNegSlope:]<(H.max()/2))[0][0]
		
			# divide by the standard deviation
			bgNoiseZ = bgNoise / bgNoise.std()
			bgChiParams = ss.chi.fit(bgNoiseZ)
			print bgChiParams
			
			# now generate values that are consistent with the histogram
			yx=range(0,H.size)/bgNoise.std()
			rvs=ss.chi.pdf(yx,bgChiParams[0],loc=bgChiParams[1],scale=bgChiParams[2])
		
			# now we can calculate the goodness of fit
			gof=np.average(np.absolute(H[halfMaxRightTail:]-rvs[halfMaxRightTail:]))
			QI2=QI1+gof
			"""

		else:
			print "QI1 can not be calculated for floating point data, skipping subject %s" % (self.sub_id)
			raise TypeError

		return QI1


pipeline_path = sys.argv[1]

csf_threshold = sys.argv[2]

gm_threshold = sys.argv[3]

wm_threshold = sys.argv[4]

subject_list = os.listdir(pipeline_path)
#subject_list = ['0025834_session_2','0025835_session_1','0025835_session_2','0025836_session_1','0025836_session_2','0025837_session_1','0025837_session_2','0025838_session_1','0025838_session_2','0025839_session_1','0025839_session_2','0025840_session_1','0025840_session_2','0025841_session_1','0025841_session_2']

'''
### CODE TO RUN ON A SINGLE CORE ###

# Function to run QC measures for all subjects and handle errors without crashing.
def run_subs(sub_list, p_path, csf_thresh, gm_thresh, wm_thresh):
	
	# Create a meta-dictionary to hold the QC dictionaries for each subject.
	qc_multi_dict = dict()

	# For each subject in the subject list...
	for sub_id in sub_list:
	
		try:
			# Create a new instance of the do_qc class, passing in the subject ID.
			subject_qc = do_qc(sub_id, p_path, csf_thresh, gm_thresh, wm_thresh)
		except:
			# If an error is encountered for this subject, print the error message and continue.
			pass
		else:
			# Add the dictionary of QC values for this subject to our meta-dictionary, with the subject ID as the key.
			# Only if the subject runs without errors.
			qc_multi_dict[sub_id] = subject_qc.qc_dict
	
	return qc_multi_dict

# Create a meta-dictionary containing info from all subjects run without multiprocessing.
qc_serial_dict = run_subs(subject_list, pipeline_path, csf_threshold, gm_threshold, wm_threshold)

# From our meta-dictionary, create a new pandas DataFrame, with meta-level keys (subjects) as rows and QC measures as columns.
qc_df = pd.DataFrame.from_dict(qc_serial_dict, orient='index')

qc_df.to_csv('/data/Projects/CoRR/preproc/qc/scripts/test.csv')

'''

### CODE TO RUN ON MULTIPLE CORES ###

# Does the same thing as run_subs, but can be passed to pool.map for multiprocessing support.
def run_subs_multiproc(sub_id, p_path, csf_thresh, gm_thresh, wm_thresh):

	try:
		# Create a new instance of the do_qc class, passing in the subject ID.
		subject_qc = do_qc(sub_id, p_path, csf_thresh, gm_thresh, wm_thresh)
	except:
		# If an error is encountered for this subject, create an empty dictonary with just the subject ID.
		sub_dict = dict()
		sub_dict['subject'] = sub_id 
		pass
	else:
		# Create a variable that allows us to return the subject's QC dictionary.
		# Only runs if no errors are encountered.
		sub_dict = subject_qc.qc_dict
		print("QC calculated for subject %s" % (sub_id))
	return sub_dict
	
# Create a "partial" function from run_subs_multiproc that sets defaults values for all arguments except for sub_id.
# This is necessary because when calling a function using pool.map, you can only specify one argument.
partial_run_subs_multiproc = partial(run_subs_multiproc, p_path=pipeline_path, csf_thresh=csf_threshold, gm_thresh=gm_threshold, wm_thresh=wm_threshold)

# Create a pool for multiprocessing
pool = Pool(48) # Use 2 cores.

# Run QC on every subject in subject_list, store the results as a list.
multiproc_store = pool.map(partial_run_subs_multiproc, subject_list)

# Close the pool and wait for the job to finish.
pool.close()
pool.join()

# Create a new dataframe from our multiproc_store, with subjects as rows and QC measures as columns.
multiproc_df = pd.DataFrame.from_records(multiproc_store).set_index('subject')

# Save the newly created DataFrame as a .csv
multiproc_df.to_csv('/data/Projects/CoRR/preproc/qc/anat_qc.csv')
