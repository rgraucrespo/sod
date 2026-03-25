#!/bin/bash
# This script runs the gcstats program of the SOD package

# Read the values of nsubsmin and nsubsmax from the file INGC
read nsubsmin nsubsmax < <(sed -n '2p' INGC)

# Copy the required OUTSOD and ENERGIES files
for ((i=nsubsmin; i<=nsubsmax; i++)); do
    cp -n ../n$(printf "%02d" $i)/OUTSOD   OUTSOD_$(printf "%02d" $i)
    cp -n ../n$(printf "%02d" $i)/ENERGIES ENERGIES_$(printf "%02d" $i)
    cp -n ../n$(printf "%02d" $i)/DATA     DATA_$(printf "%02d" $i) 2>/dev/null || true
    cp -n ../n$(printf "%02d" $i)/SPECTRA  SPECTRA_$(printf "%02d" $i) 2>/dev/null || true
    cp -n ../n$(printf "%02d" $i)/XSPEC    . 2>/dev/null || true
done

gcstatsod

rm -f OUTSOD* ENERGIES* SPEC* DATA*
