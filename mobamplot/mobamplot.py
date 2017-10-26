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
import argparse
import inspect
import os
import subprocess
import warnings
import sys
import h5py
import numpy as np
import scipy.stats as stats
import statistics
import math

import plotly
import plotly.plotly as py
import plotly.tools as tls
import plotly.graph_objs as go
import plotly.plotly as py

from plastid.util.io.openers import get_short_name
from plastid.util.io.filters import NameDateWriter
from plastid.util.scriptlib.help_formatters import format_module_docstring


# ignore and do not print any occurences of matching warnings
warnings.simplefilter("ignore")
# A stream to which stderr-like info can be written
printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))

# call getCoutVectorData to create hdf5 file with counts
def callFunc(bedfile, bamList, bamIDs, outfolder):
	printer.write("Calling getcountvectordata.py")

	# create list of arguments to feed into fxn as if calling from the terminal
	args = ['python', '/home/monika/myBamPlot/bamPlotCode/getCountVectorData.py', bedfile, '--bamfiles'] + bamList + ['--outfolder', outfolder, '--bamIDs'] + bamIDs
	
	# call function
	process = subprocess.Popen(args, stdout=subprocess.PIPE, universal_newlines=True)
	# print stdout messages of fxn
	print(process.communicate()[0].encode('utf-8').decode('unicode_escape'))

# open dsets in file, 38
def plotDataHist(hdf5, filetype, outfolder, viewplot, binsize):
	# open hdf5 file
	f = h5py.File(hdf5, "r")

	# extract data from file
	for ID, dset in f.items():

		chrom = dset.attrs["chrom"].decode("utf-8")		# chromosome
		chromStart = dset.attrs["chromStart"]			# start coordinate
		chromEnd = dset.attrs["chromEnd"]				# end coordinate
		x_data = np.arange(chromStart, chromEnd)		# position np.array (end exclusive)

		bamIDs = list(dset.dtype.names)					# list of bamIDs

		subplots = list(range(1, len(dset.dtype)+1))	# list of subplot numbers

		# create subplot
		fig = tls.make_subplots(rows=len(subplots), subplot_titles = bamIDs,
								vertical_spacing=0.05)

		# Update 'data' key in fig with a Histogram object for every BAM to respective subplot
		fig['data'] = go.Data([make_trace_hist(x_data, dset[BAM], BAM, binsize, index) for index, BAM in enumerate(bamIDs, start=1)])

		# specify layout
		fig = setLayout(fig, chrom, chromStart, chromEnd, subplots, ID)
		fig = addSlider(fig)

		# files saved to specified output folder and named by ROI
		if filetype == None:
			plot_url = plotly.offline.plot(fig, filename=str(outfolder) + "/" + str(chrom) + "_.hist", 
										   auto_open=viewplot)
		elif (filetype == 'svg') or (filetype == 'webp'):
			# saves file types [png, jpeg, svg, webp] AND HTML, but must open plot in web browser
			plot_url = plotly.offline.plot(fig, filename=str(outfolder) + "/" + str(chrom) + "_.hist",
											image=filetype, image_filename=str(chrom) + "_.hist",
											auto_open=viewplot)
		else:
			# automatically saves files of type [png, jpeg] without opening browser
			# No HTML file saved
			py.image.save_as(fig, filename=str(outfolder) + "/" + str(chrom) + "_.hist" + filetype)

	f.close()
	printer.write("Figure completed.")

# plot data function
def plotDataLine(hdf5, filetype, outfolder, viewplot, binsize):

	# open hdf5 file
	f = h5py.File(hdf5, "r")

	# extract data from file
	for ID, dset in f.items():

		chrom = dset.attrs["chrom"].decode("utf-8")		# chromosome
		chromStart = dset.attrs["chromStart"]			# start coordinate
		chromEnd = dset.attrs["chromEnd"]				# end coordinate
		x_data = np.arange(chromStart, chromEnd)		# position np.array (end exclusive)

		### BINNING
		# create array of values of nbin spaces from chromstart to closest evenly divisible value to chromEnd
		x_data_binned = np.arange(chromStart, chromEnd, binsize)
		x_data_binned = np.append(x_data_binned, chromEnd-1)
		# list x_data_binned centers for plotting
		x_data_binned_mean = [statistics.mean([x_data_binned[i], x_data_binned[i+1]]) for i in range(0, len(x_data_binned)-1)]

		bamIDs = list(dset.dtype.names)				# list of bamIDs

		subplots = list(range(1, len(dset.dtype)+1))	# list of subplot numbers

		# create subplot
		fig = tls.make_subplots(rows=len(subplots), vertical_spacing=0.05) #subplot_titles = bamIDs,

		# Update 'data' key in fig with a Histogram object for every BAM to respective subplot
		fig['data'] = go.Data([make_trace_line(x_data, x_data_binned, x_data_binned_mean, dset[BAM], BAM, index) for index, BAM in enumerate(bamIDs, start=1)])

		# specify layout
		fig = setLayout(fig, chrom, chromStart, chromEnd, subplots, ID)

		# files saved to specified output folder and named by ROI
		if filetype == None:
			plot_url = plotly.offline.plot(fig, filename=str(outfolder) + "/" + str(chrom) + "_.line.bin" + str(binsize), 
										   auto_open=viewplot)
		elif (filetype == 'svg') or (filetype == 'webp') or (filetype == 'png') or (filetype == 'jpeg'):
			# saves file types [png, jpeg, svg, webp] AND HTML, but must open plot in web browser
			plot_url = plotly.offline.plot(fig, filename=str(outfolder) + "/" + str(chrom) + "_.line.bin" + str(binsize),
											image=filetype, image_filename=str(chrom) + "_.line.bin" + str(binsize),
											auto_open=viewplot, image_height=1200, image_width=800)
		else:
			# automatically saves files of type [png, jpeg] without opening browser
			# No HTML file saved
			py.image.save_as(fig, filename=str(outfolder) + "/" + str(chrom) + "_.line.bin" + str(binsize) + filetype)

	f.close()
	printer.write("Figure completed.")

		# finds the mean of the binned counts when ROI split into bin_size # of bins
		# gets mean values for binned x_data, this is point that is used for plotting
		#bin_means, bin_edges, binnumber = stats.binned_statistic(x_data, counts[BAM], statistic='mean', bins=binsize)
		#x_data_binned_mean, xbin_edges, xbinnumber = stats.binned_statistic(x_data, x_data, statistic='mean', bins=binsize)

# histogram trace-generating function
def make_trace_hist(x_data, counts, BAM, binsize, index):
	trace = go.Histogram(
				x=x_data,					# x-coords are same for graphs in a ROI
				y=counts,					# y-coords are counts for specified BAM
				histfunc="avg",				# specifies binning function, values from avg in bin
				autobinx=False,				# disables x axis automatic binning
				xbins=go.XBins(
					start=x_data[0],		# set start value for bins
					end=x_data[-1],			# set end value for bins (1 more than counts to include final point in calcs)
					size=binsize),			# set step value in-between each bin
				name=BAM,					# set trace name to respective BAM
				xaxis="x{}".format(index),	# references the respective subplot
				yaxis="y{}".format(index),	# references the respective subplot
				marker=dict(
					color="0000ff"
					)
			)
	return trace

# line trace-generating function
def make_trace_line(x_data, x_data_binned, x_data_binned_mean, counts, BAM, index):

	# countAvg = binned count values, bin_edges = returns x_data binned edges, binnumber = ndarray of which bin each count value corresponds to
	countAvg, bin_edges, binnumber = stats.binned_statistic(x_data, counts, statistic='mean', bins=x_data_binned)

	trace = go.Scatter(
				x=x_data_binned_mean,		# x-coords are same for graphs in a ROI
				y=countAvg,					# y-coords are avg counts/bin for specified BAM
				fill="tozeroy",				# color fill to y=0
				fillcolor="0000ff",
				mode="line",				# line plot
				opacity=1,
				line=dict(
					color="0000ff",
					shape='linear'),		# lines drawn using spline interpolation
					#smoothing=0),			# amount of smoothing
				name=BAM,
				xaxis="x{}".format(index),	# references the respective subplot
				yaxis="y{}".format(index)	# references the respective subplot
			)
	return trace

# set Layout features of graph
def setLayout(fig, chrom, chromStart, chromEnd, subplots, ROI_name):
	# title specs
	titlefont = dict(size=18)

	# Define dictionary of subplot xaxis layout
	xLay = dict(
		range=[chromStart-300,chromEnd+300],			# graph range larger than ROI
		showline=True,									# draw line bounding axes
		tickmode="array",								# placement of ticks set w/ tickvals
		tickvals=[chromStart,chromEnd],					# xaxis values of text
		ticktext=[str(chrom) + ": " + str(chromStart), 	# values of text, coords at ROI start & end
				  str(chrom) + ": " + str(chromEnd)],
		ticks="outside",								# draw ticks outside axes
		zeroline=True,									# draw line along 0 value of axis
		zerolinewidth=2,								# set zero line width
		mirror=True										# axis lines mirrored (also wide)
	)
	# max y-axis value of all subplots
	countsMax = max([max(fig.data.get_data()[i]['y']) for i in range(0, subplots[-1])])

	# round countsMax up to nearest 0.5 for yrange max
	y_range_max = 0.5 * math.ceil(2 * countsMax)

	# look @ BAM liquidator to see how it's getting counts
	# difference in counts (read extension?)
	# ceiling function, entire array, find max and then wrap all into ceiling function to closest highest integer to prevent fractional values as max

	# Define dictionary of subplot yaxis layout
	yLay = dict(
		title="Reads",				# y-axis title "Reads"
		range=[0, y_range_max],		# graph range
		tick0=0,					# placement of first tick		
		dtick=(y_range_max/20),					# step size between ticks
		ticks="outside",			# draw ticks outside axes
		showgrid=False,				# omit grid
		showline=True,				# show line bounding axes
		zeroline=True,				# draw line along y=o
		zerolinewidth=2,			# set zero line width
		mirror=True					# axis lines mirrored (also wide)
	)

	legend = dict(
		bordercolor="rgba(0,0,0,1)",
		#x=0.5,
		#y=1,
		orientation="h"	
	)

	for splt in subplots:
		# Update x & y axis layout of all subplots
		fig["layout"].update({"xaxis{}".format(splt): xLay})
		fig["layout"].update({"yaxis{}".format(splt): yLay})
		fig["layout"].update(legend=legend)

		# create scale bar for each subplot
		bar_labels = go.Annotation(text="1kb", 								# text above scalebar
					   x=chromEnd-500, y=y_range_max*(.923), 				# x at center of bar, y in upper right corner
					   xref="x{}".format(splt), yref="y{}".format(splt), 	# set subplot
					   showarrow=False, 									# no arrow
					   xanchor="center", yanchor="bottom") 					# text centered at x, bottom of box at y

		bars = dict(type="rect", layer="above", 							# rectangle, shape drawn above trace
				 xref="x{}".format(splt), yref="y{}".format(splt), 			# set subplot
				 x0=chromEnd-1000, x1=chromEnd, 							# 1000 bp wide
				 y0=y_range_max*(0.903), y1=y_range_max*(0.923), 			# height
				 fillcolor="rgba(0,0,0,1)") 								# fill black

		# update layout with scale bar
		fig["layout"]["shapes"].append(bars)
		fig["layout"]["annotations"].append(bar_labels)

	# add graph title
	fig["layout"].update(title="Gene: " + ROI_name, titlefont=titlefont, showlegend=False, autosize=True)
	return fig

# add slider option, only for hist
def addSlider(fig):

	# create slider
	steps = []
	for i in np.arange(10, 200, 10): 	# from 10 to 200 in step size of 10
		step = dict(
			method="restyle",			# "restyle": modifying data attribute
			args = ['xbins.size', i],	# spec changing attribute, bin values
			label = str(i)				# label step as size
		)
		steps.append(step)

	sliders = [dict(
		currentvalue = {"prefix": "Bin Size: "},	# prefix of current slider value label
		pad = {"t": 15},								# 50 px along top of slider
		steps = steps
	)]

	# add slider to layout
	fig['layout'].update(sliders=sliders)
	return fig

# determines color of plot
def prettyColors(a):
	colorDict = dict(
			blue="rgba(0,0,255,{})".format(a),
			red="rgba(255,0,0,{})".format(a),
			yellow="rgba(255,255,0,{})".format(a),
			green="rgba(118,92,46,{})".format(a),
			purple="rgba(278,68,47,{})".format(a),
			aqua="rgba(181,100,47,{})".format(a),			
			orange="rgba(26,100,47,{})".format(a),
			pink="rgba(312,100,57,{})".format(a)
		)

	colorList = []

# meta option for plot, req group assignment which should be plotted within the same subplot as a meta
# have a separate plotting function for plotting metas
# current meta plots (feed list of BAM, for each BAM string (same length) that's associated

# bam option to feed folder, call function from general pipeline script with designated data table

# profile for bottle necks
# python docs module (read the docs)

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
	parser = argparse.ArgumentParser(description=format_module_docstring(__doc__), 
									 formatter_class=argparse.RawDescriptionHelpFormatter)

	# required files
	parser.add_argument('--bedfile', type=str, help="Input BED file.")
	parser.add_argument('--hdf5', type=str, default=None,
						 help="Input hdf5 file with counts.")

	# required files
	bamgroup = parser.add_mutually_exclusive_group(required=False)
	bamgroup.add_argument('--bamfolder', help="Input folder containing BAM files.")
	bamgroup.add_argument('--bamfiles', type=str, nargs='+', help="Input BAM files.")

	# folder to save files/figures
	parser.add_argument('--outfolder', type=str, help="Folder to save the HTML plots. \
														If none is stated files will be saved to current working directory.")
	# optinal args
	parser.add_argument('--plot', action='store_true', default=False, help="Plot data to figure.")
	parser.add_argument('--viewplot', action='store_true', default=False, help="Open web browser to view plot.")
	parser.add_argument('--binsize', type=int, default=200, help="Bin size.")
	parser.add_argument("--format", default=None, type=str, choices=[None, 'png', 'svg', 'jpeg', 'webp'],
						help="Format for static image file to be saved as.")
	parser.add_argument("--line", default=False, action='store_true', help="Plot as a filled line graph.")
	parser.add_argument("--bamIDs", default=None, nargs='+',
						help="Identification for each BAM that will displayed for its respective plot.")
	parser.add_argument("--callFunc", default=False, action='store_true', help="Calls getCountVectorData.py.")
	parser.add_argument("--metaplot", default=False, action='store_true', help="Create a metaplot.")

	args = parser.parse_args()

	# parse args
	bedfile = args.bedfile
	hdf5 = args.hdf5
	outfolder = args.outfolder

	viewplot = args.viewplot
	binsize = args.binsize
	filetype = args.format

	# creates list of BAMs from folder or given files
	if (args.bamfolder != None):
		bamList = [file for file in glob.glob(os.path.join(args.bamfolder, "*.bam*"))]
	else:
		bamList = args.bamfiles

	#if (args.bamIDs != None) and (len(args.bamIDs) != len(bamList)):
	#	print("InputError: Number of BAM_IDs MUST match number of BAMS. \n\
	 #   number of BAM_IDs: {0}  number of BAMs: {1}".format(len(args.bamIDs), len(bamList)))
	#	exit()

	if (bamList != None):
		if (args.bamIDs == None):
		# List with short BAM names (without filepath or extention)
			bamIDs = [os.path.splitext(os.path.basename(bamfile))[0] for bamfile in bamList]
		else:
			bamIDs = args.bamIDs

	# if output folder not specified, use cwd
	# if output folder doesn't exist, create it
	# print output folder with path
	if outfolder == None:
		outfolder = os.getcwd()
		printer.write("No output folder specified.")
	else:
		if not os.path.isdir(outfolder):
			os.mkdir(outfolder)
			printer.write("Output folder does not exist. Created folder.")
	printer.write("Saving plot output files to: {}".format(os.path.realpath(outfolder)))	

	if args.callFunc == True:
		hdf5 = callFunc(bedfile, bamList, bamIDs, outfolder)

	if args.line == False:
		plotDataHist(hdf5, filetype, outfolder, viewplot, binsize)
	else:
		plotDataLine(hdf5, filetype, outfolder, viewplot, binsize)
	return

if __name__ == "__main__":
	main()
