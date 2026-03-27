#!/bin/bash
# This script runs the spbe (simple pair-based extrapolation) programe of the SOD package.
# It requires the following input files: EQMATRIX, OUTSOD, OUTSOD1, OUTSOD2, ENERGIES0, ENERGIES1, ENERGIES2
# The following commands might be useful to collect the input files - comment out if you already have them.
cp ../../EQMATRIX ./
cp ../OUTSOD ./
cp ../../n01/OUTSOD ./OUTSOD1
cp ../../n02/OUTSOD ./OUTSOD2
cp ../../n00/ENERGIES ./ENERGIES0
cp ../../n01/ENERGIES ./ENERGIES1
cp ../../n02/ENERGIES ./ENERGIES2

spbesod

rm -f EQMATRIX OUTSOD* ENERGIES? 
