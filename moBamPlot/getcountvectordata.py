#!/usr/bin/env python
"""Fetch vectors of :term:`counts` at each nucleotide position in one or more
regions of interest (ROIs) and constructs a list of tuples for each ROI with
their respective BAM-associated count_vectors.


Output files
------------
HTML files save in a user-specified output folder or current working directoy as individual files 
named after the ROI's chromosome and the binsize to which it corresponds.
HTML files do not automatically open in current web browswer or open in new webpage with one ROI per tab.
"""

### python getCountVectorData.py ../skeratinocyte_genes.bed --bamfiles /media/monika/My\ Passport/bams/Yng_H3K27ac_1-1.hg19.bwt2.sorted.bam

import argparse
import inspect
import os
import warnings
import sys
import h5py
import glob
import numpy as np

from plastid.genomics.genome_array import BAMGenomeArray
from plastid.util.io.openers import get_short_name
from plastid.util.io.filters import NameDateWriter
from plastid.util.scriptlib.help_formatters import format_module_docstring
from plastid.genomics.roitools import SegmentChain
from plastid.readers.bed import BED_Reader

# ignore and do not print any occurences of matching warnings
warnings.simplefilter("ignore")

# A stream to which stderr-like info can be written
# when printed in terminal, includes the script name, time & date
printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))

# create list of tuples [(SegmentChain, counts)]
def getCountVectorData(bedfile, bamList, bamIDs):
	
	# read BED file to iterator of SegmentChain objects
	bed_segmentChains = list(BED_Reader(open(bedfile)))

	# create list of tuples (SegmentChain, [count_vectors/BAM])
	countVectorData = []

	# iterate through ROIs
	for index,segment in enumerate(bed_segmentChains):

		# prints in terminal progress every 100 ROIs processed
		if index % 100 == 0:
			#printer.write("Processed %s regions of interest" % index)
			print("Processed %s regions of interest" % index)

		# create generator of BAMGenomeArray objects
		print(bamList)
		BAM_al = (BAMGenomeArray(bamfile) for bamfile in bamList)

		# list of count vectors per BAM
		count_vec = [segment.get_counts(al) for al in BAM_al]

		# create recorded array, columns are count vectors
		# MUST call by BAM name, NOT ind
		count_array = np.core.records.fromarrays(count_vec, names=bamIDs)

		# add tuple to countVectorData list
		countVectorData.append((segment, count_array))

	#printer.write("Completed processing data.")
	print("Completed processing data.")
	return countVectorData

# save data to hdf5 file
def saveHDF5(bedfile, bamList, bamIDs, outfolder):

	countVectorData = getCountVectorData(bedfile, bamList, bamIDs)

	# create h5py file in output folder
	file = h5py.File(outfolder + "/" + os.path.splitext(os.path.basename(bedfile))[0] + ".hdf5", "w")

	for segment, counts in countVectorData:

		# create dataset
		dset = file.create_dataset(str(segment.get_name()), data=counts)

		# add segment chain information to dataset as attributes
		dset.attrs["ID"] = np.string_(segment.get_name())
		dset.attrs["chrom"] = np.string_(segment.chrom)
		dset.attrs["chromStart"] = segment.segments[0].start
		dset.attrs["chromEnd"] = segment.segments[0].end

	# record absolute filepath, then close file
	filepath = os.path.abspath(file.filename)
	file.close()
	#printer.write("Data successfully saved to .hdf5 file.")
	print("Data successfully saved to .hdf5 file.")
	return filepath

def main(args=sys.argv[1:]):
	"""Command-line program
	
	Parameters
	----------
	argv : list, optional
		A list of command-line arguments, which will be processed
		as if the script were called from the command line if
		:func:`main` is called directly.

		Default: `sys.argv[1:]`. The command-line arguments, if the script is
		invoked from the command line
	"""
	# create ArgumentParser object to hold all info to parse the command line into Python data types
	# description from docstring in exact format, format_module_docstring pretty prints docstrings
	parser = argparse.ArgumentParser(description="Makes dat data tabullllll", 
									 formatter_class=argparse.RawDescriptionHelpFormatter)

	# required files
	parser.add_argument('bedfile', type=str, help="Input BED file.")
	group = parser.add_mutually_exclusive_group(required=True)
	group.add_argument('--bamfolder', help="Input folder containing BAM files.")
	group.add_argument('--bamfiles', type=str, nargs='+', help="Input BAM files.")

	# folder to save files/figures
	parser.add_argument('--outfolder', type=str, help="Folder to save the HTML plots. \
														If none is stated files will be saved to current working directory.")
	# optinal args
	parser.add_argument("--bamIDs", default=None, nargs='+',
						help="Identification for each BAM.")

	args = parser.parse_args()
	
	# parse args
	bedfile = args.bedfile
	outfolder = args.outfolder

	# creates list of BAMs from folder or given files
	if (args.bamfolder != None):
		bamList = [file for file in glob.glob(os.path.join(args.bamfolder, "*.bam"))]
	else:
		bamList = args.bamfiles

	if (args.bamIDs != None) and (len(args.bamIDs) != len(bamList)):
		print("InputError: Number of BAM_IDs MUST match number of BAMS. \n\
		number of BAM_IDs: {0}  number of BAMs: {1}".format(len(args.bamIDs),len(bamList)))
		exit()

	if (args.bamIDs == None):
		# List with short BAM names (without filepath or extention)
		bamIDs = [os.path.splitext(os.path.basename(bam_file))[0] for bam_file in bamList]
	else:
		bamIDs = args.bamIDs

	# if output folder not specified, use cwd
	# if output folder doesn't exist, create it
	# print output folder with path
	if outfolder == None:
		outfolder = os.getcwd()
		print("No output folder specified.")
	else:
		if not os.path.isdir(outfolder):
			os.mkdir(out_older)
			print("Output folder does not exist. Created folder.")
	#printer.write("Saving counts data files to: {}".format(os.path.realpath(outfolder)))
	print("Saving counts data files to: {}".format(os.path.realpath(outfolder)))
	# run saveHDF5
	saveHDF5(bedfile, bamList, bamIDs, outfolder)


if __name__ == "__main__":
	main()
