"""A tool to process and graph bam read density at specific loci.

"""
import os

from setuptools import setup, find_packages

setup(
  name='moBamPlot',
	version='0.1dev',
	description='A Python program to plot bams.',
	long_description = open('README.txt').read(),

	# LinLab at BCM
	url='PUT GITHUB URL HERE',
	download_url='PUT DOWNLOAD URL TO GITHUB HERE',

	author='Monika Perez',
	author_email='mp22@bcm.edu',

	packages=find_packages(),

	install_requires=[],
	extras_require={},

	scripts=['scripts/getCountVectorData.py','scripts/moBamPlot.py']
)
