"""A tool to process and graph bam read density at specific loci.
"""
import os
from setuptools import setup, find_packages
from setuptools.command.install import install
from subprocess import check_call

setup(
  name='mobamplot',
	version='0.1dev0',
	description='A Python program to plot bams.',
	long_description = open('README.md').read(),

	# LinLab at BCM
	url='https://github.com/linlabcode/moBamPlot',

	author='Monika Perez',
	author_email='mp22@bcm.edu',

	packages=find_packages(),

	install_requires=[
		"numpy >= 1.9.0",
		"Pysam >= 0.8.4",
		"Cython >= 0.22.0",
		"h5py == 2.7.0",
		"scipy == 0.19.0",
		"plotly == 2.0.15"
	],
	extras_require={},

	scripts=['moBamPlot/getCountVectorData.py','moBamPlot/moBamPlot.py']
)
