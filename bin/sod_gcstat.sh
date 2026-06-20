#!/bin/bash
# This script runs the gcstats program of the SOD package
# Run from an x????/ folder at the same level as the nXX/ folders

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="${SCRIPT_DIR}:${PATH}"
. "${SCRIPT_DIR}/sod_common.sh"

SODPROJECT="$(sod_require_project_root "$PWD")" || exit 1
XDIR="$(sod_find_ancestor_with_file "$SODPROJECT" "$PWD" INGC)" || {
    echo "Error: INGC not found in the current folder or its parents below SODPROJECT/."
    exit 1
}
cd "$XDIR" || exit 1

# Read the values of nsubsmin and nsubsmax from the file INGC
read nsubsmin nsubsmax < <(sed -n '2p' INGC)

# Copy the required ENSEMBLE and ENERGIES files
for ((i=nsubsmin; i<=nsubsmax; i++)); do
    ndir="$(sod_level_dir_by_number "$SODPROJECT" "$i")" || {
        echo "Error: level directory for nsubs=$i not found in SODPROJECT/."
        exit 1
    }
    tag="${ndir#n}"
    cp -n "$SODPROJECT/$ndir/ENSEMBLE" ENSEMBLE_${tag}
    cp -n "$SODPROJECT/$ndir/ENERGIES" ENERGIES_${tag}
    cp -n "$SODPROJECT/$ndir/DATA"     DATA_${tag}   2>/dev/null || true
    cp -n "$SODPROJECT/$ndir/SPECTRA"  SPECTRA_${tag} 2>/dev/null || true
    cp -n "$SODPROJECT/$ndir/XSPEC"    . 2>/dev/null || true
done

gcstatsod

rm -f ENSEMBLE* ENERGIES* SPEC* DATA*
