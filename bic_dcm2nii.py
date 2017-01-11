"""
BIC DICOM to NIfTI Conversion Script - Written by Dan Lurie (dan.lurie@berkeley.edu)

Designed for use on data from the UC Berkeley Brain Imaging Center.

Requires dcm2nii (http://www.mccauslandcenter.sc.edu/mricro/mricron/dcm2nii.html).

This script will convert all anatomical (MPRAGE) and functional (EPI) scans for all subjects in a directory.

Usage:
======

python bic_dcm2nii dicom_dir nifti_dir


Notes:
====== 
   
dicom_dir: The path to a directory containing your subject folders.
           Each subject folder should contain the following:
	   - At least one anatomical scan folder.
	   - At least one functional scan folder.

nifti_dir: A directory where you would like NIfTI files written.
	   If this directory does not exist, one will be created.

anat_name: A string common to all anatomical scan folder names.
           Example: Given folders MPRAGE_1 & MPRAGE_2, anat_name = MPRAGE

func_name: A string common to all functional scan folder names.
           Example: Given folders Rest1 & Rest2, func_name = Rest
"""

import os
import sys
import subprocess

dicom_dir = sys.argv[1]
nifti_dir = sys.argv[2]
anat_name = sys.argv[3]
func_name = sys.argv[4]

subject_list = os.listdir(dicom_dir)

for subject in subject_list:
    # Get the path to the raw data for this subject.
    subject_raw_path = '/'.join([dicom_dir, subject])
    # Set the directory where NIfTIs will be written for this subject.
    subject_nifti_path = '/'.join([nifti_dir, subject])
    # Set paths for where we will store anatomical and functional NIfTIs.
    func_path = '/'.join([subject_nifti_path, 'func'])
    anat_path = '/'.join([subject_nifti_path, 'anat'])
    # Create directories for anatomical and functional NIfTIs.
    if os.path.isdir(func_path) == False:
        os.makedirs(func_path)
    else:
        pass
    if os.path.isdir(anat_path) == False:
        os.makedirs(anat_path)
    else:
        pass
    # Get a list of all scans we have for this subject.
    scan_list = os.listdir(subject_raw_path)
    # Get a list of all MPRAGE scans.
    mprage_list = filter(lambda x: anat_name in x, scan_list)
    # Get a list of all EPI scans.
    epi_list = filter(lambda x: func_name in x, scan_list)
    # Convert MPRAGE scans
    for mprage_scan in mprage_list:
        # Get the path of the DICOMs for this scan.
        mprage_raw_path = '/'.join([subject_raw_path, mprage_scan])
        # Create a folder in ../anat for this scan.
        mprage_nifti_dir = '/'.join([anat_path, mprage_scan])
        if os.path.isdir(mprage_nifti_dir) == False:
            os.makedirs(mprage_nifti_dir)
            # Call dcm2nii to do the actual conversion.
            subprocess.check_call(['dcm2nii', '-o', mprage_nifti_dir, mprage_raw_path])
        else:
            pass
    # Convert EPI scans
    for epi_scan in epi_list:
        # Get the path of the DICOMs for this scan.
        epi_raw_path = '/'.join([subject_raw_path, epi_scan])
        # Create a folder in ../anat for this scan.
        epi_nifti_dir = '/'.join([func_path, epi_scan])
        if os.path.isdir(epi_nifti_dir) == False:
            os.makedirs(epi_nifti_dir)
            # Call dcm2nii to do the actual conversion.
            subprocess.check_call(['dcm2nii', '-o', epi_nifti_dir, epi_raw_path])
        else:
            pass	
