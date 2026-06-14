#!/bin/bash
# sod_mcstat.sh — Thermodynamic integration over MC temperatures.
#
# Run from nXX/ (where the MCT_*K/ directories reside).
# Reads E_ave from each MCT_TTTK/PMEx/ENSEMBLE + MCT_TTTK/PMEx/ENERGIES and
# performs Gibbs-Helmholtz thermodynamic integration with the T=∞ reference
# S(∞) = kB ln(C(npos,lev)).  The PMEx variant is read from ../pme.model when
# present, otherwise it defaults to PMEh (matching mcsod's default).
#
# Output: thermodynamics.dat  (same format as statsod/gcstatsod)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="${SCRIPT_DIR}:${PATH}"

if [ ! -f ../TEMPERATURES ]; then
  echo "Error: required file '../TEMPERATURES' not found."
  echo "Run sod_mcstat.sh from nXX/ (where the MCT_*K/ directories are)."
  exit 1
fi

if ! ls -d MCT_*K/ 2>/dev/null | grep -q .; then
  echo "Error: no MCT_*K/ directories found. Run sod_mc.sh first."
  exit 1
fi

mcstatsod
