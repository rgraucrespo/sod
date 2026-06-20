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
. "${SCRIPT_DIR}/sod_common.sh"

SODPROJECT="$(sod_require_project_root "$PWD")" || exit 1
LEVEL_NAME="$(sod_find_enclosing_level_name "$SODPROJECT" "$PWD" || true)"
if [ -z "$LEVEL_NAME" ]; then
  echo "Error: run sod_mcstat.sh from an nXX/ folder or one of its subfolders."
  exit 1
fi
cd "$SODPROJECT/$LEVEL_NAME" || exit 1

if [ ! -f ../TEMPERATURES ]; then
  echo "Error: required file '../TEMPERATURES' not found."
  echo "Run sod_mcstat.sh from nXX/ or one of its subfolders."
  exit 1
fi

if ! ls -d MCT_*K/ 2>/dev/null | grep -q .; then
  echo "Error: no MCT_*K/ directories found. Run sod_mc.sh first."
  exit 1
fi

mcstatsod
