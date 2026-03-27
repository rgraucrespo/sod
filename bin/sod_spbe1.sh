#!/bin/bash
# This script runs the spbe code (simple pair-based extrapolation) programe of the SOD package, starting from x=1.
# It requires the following input files: EQMATRIX, OUTSOD, OUTSOD1, OUTSOD2, ENERGIES0, ENERGIES1, ENERGIES2
# The following commands might be useful to collect the input files - comment out if you already have them.
nsites=`head -1 ../OUTSOD | awk '{print $4}'`
nsitesm1=$(($nsites - 1))
nsitesm2=$(($nsites - 2))

# make sure that nsites has the right number of digits
if (( nsites < 10 )); then
    nsites=$(printf "%02d" "$nsites")
elif (( nsitesm1 < 100 )); then
    nsites=$(printf "%03d" "$nsites")
else
    nsites=$(printf "%04d" "$nsites")
fi

# make sure that nsitesm1 has the right number of digits
if (( nsitesm1 < 10 )); then
    nsitesm1=$(printf "%02d" "$nsitesm1")
elif (( nsitesm1 < 100 )); then
    nsitesm1=$(printf "%03d" "$nsitesm1")
else
    nsitesm1=$(printf "%04d" "$nsitesm1")
fi

# make sure that nsitesm2 has the right number of digits
if (( nsitesm2 < 10 )); then
    nsitesm2=$(printf "%02d" "$nsitesm2")
elif (( nsitesm2 < 100 )); then
    nsitesm2=$(printf "%03d" "$nsitesm2")
else
    nsitesm2=$(printf "%04d" "$nsitesm2")
fi

cp ../../EQMATRIX ./
cp ../OUTSOD ./OUTSOD_original
invertOUTSOD
mv OUTSOD_inverted OUTSOD
cp ../../n$nsitesm1/OUTSOD ./OUTSOD_original
invertOUTSOD
mv OUTSOD_inverted OUTSOD1 
cp ../../n$nsitesm2/OUTSOD ./OUTSOD_original
invertOUTSOD
mv OUTSOD_inverted OUTSOD2 


cp ../../n$nsites/ENERGIES   ./ENERGIES0
cp ../../n$nsitesm1/ENERGIES ./ENERGIES1
cp ../../n$nsitesm2/ENERGIES ./ENERGIES2
rm ./OUTSOD_original
spbesod

rm -f EQMATRIX OUTSOD* ENERGIES?

