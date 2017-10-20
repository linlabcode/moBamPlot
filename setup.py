"""A tool to process and graph bam read density at specific loci.
"""
import os

from setuptools import setup, find_packages


setup(
  name='mobamplot',
	version='0.1dev',
	description='A Python program to plot bams.',
	long_description = open('README.md').read(),

	# LinLab at BCM
	url='https://github.com/linlabcode/moBamPlot',

	author='Monika Perez',
	author_email='mp22@bcm.edu',

	packages=find_packages(),
	
	setup_requires=[
		"pysam == 0.8.4",
		"cython == 0.22.0",
		"numpy == 1.12.1"
	],

	install_requires=[
		"pysam == 0.8.4",
		"cython == 0.22.0",
		"numpy == 1.12.1",
		"h5py == 2.7.0",
		"scipy == 0.19.0",
		"plotly == 2.0.15",
		"plastid == 0.4.8"
	],
	extras_require={},

	scripts=['moBamPlot/getCountVectorData.py','moBamPlot/moBamPlot.py']
)
