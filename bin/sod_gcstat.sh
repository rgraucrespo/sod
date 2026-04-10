#!/bin/bash
# This script runs the gcstats program of the SOD package
# Run from an x????/ folder at the same level as the nXX/ folders

# Read the values of nsubsmin and nsubsmax from the file INGC
read nsubsmin nsubsmax < <(sed -n '2p' INGC)

# Determine zero-padding from existing nXX folders
sample=$(ls -d ../n*/ 2>/dev/null | head -1)
nlen=$(basename "$sample" | tr -cd '0-9' | wc -c | tr -d ' ')
[ -z "$nlen" ] || [ "$nlen" -lt 2 ] && nlen=2
fmt="%0${nlen}d"

# Copy the required OUTSOD and ENERGIES files
for ((i=nsubsmin; i<=nsubsmax; i++)); do
    tag=$(printf "$fmt" $i)
    cp -n ../n${tag}/OUTSOD   OUTSOD_${tag}
    cp -n ../n${tag}/ENERGIES ENERGIES_${tag}
    cp -n ../n${tag}/DATA     DATA_${tag}   2>/dev/null || true
    cp -n ../n${tag}/SPECTRA  SPECTRA_${tag} 2>/dev/null || true
    cp -n ../n${tag}/XSPEC    . 2>/dev/null || true
done

gcstatsod

rm -f OUTSOD* ENERGIES* SPEC* DATA*
