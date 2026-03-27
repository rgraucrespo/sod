#!/bin/bash
# This script runs the fortran code peaks2spec, which converts a list of peak positions to a spectrum (intensity vs x) by broadening with Gaussians.
#
# It requires two input files (with fixed names):
# PEAKS contains a list of peaks in each line; each line represents a different configuration (a different spectrum is generated for each configuration)
# INP2S contains other info needed to generate the spectra, e.g. xmin, xmax, broadening(sigma), etc. See an example for the format.
#
# It generates two output files:
# SPECTRA, where each line contains the generated intensity values in the x grid, for each configuration.
# XSPEC, which contains the list of x values at which the intensities are given.

peaks2spec

