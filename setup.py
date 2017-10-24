#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A tool to process and graph bam read density at specific loci."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

setup(
  name='mobamplot',
	version='0.1d.0',
	description='A Python program to plot bams.',
	long_description = readme,

	# LinLab at BCM
	url='https://github.com/linlabcode/mobamplot',

	author='Monika Perez',
	author_email='mp22@bcm.edu',

	packages=find_packages(include=['mobamplot']),
	
	include_package_data=True,
	install_requires=[
		"numpy >= 1.9.0",
		"Pysam >= 0.8.4",
		"Cython >= 0.22.0",
		"h5py == 2.7.0",
		"scipy == 0.19.0",
		"plotly == 2.0.15"
	],

	scripts=['moBamPlot/getcountvectordata.py','moBamPlot/mobamplot.py']
)
