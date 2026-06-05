#!/bin/bash
# Run from MAIN/ to score GQS at all temperatures for all nXX/ compositions, or from a single nXX/ folder.
# Reads ENERGIES from nXX/ folder, TEMPERATURES from nXX/ (priority) or MAIN/.
# Automatically runs sqssod first if OUTSQS is missing.

SOD_BIN=$(cd "$(dirname "$0")" && pwd)

check_main() {
  if [ ! -f EQMATRIX ]; then
    echo "Error: EQMATRIX not found in $(pwd). Run combsod first."
    exit 1
  fi
}

run_for() {
  local nxx="$1"
  if [ ! -f "${nxx}/ENSEMBLE" ]; then
    echo " > Skipping ${nxx}: ENSEMBLE not found."
    return
  fi
  if [ ! -f "${nxx}/ENERGIES" ]; then
    echo " > Skipping ${nxx}: ENERGIES not found."
    return
  fi
  if [ ! -f "${nxx}/OUTSQS" ]; then
    echo " > OUTSQS not found for ${nxx}. Running sqssod first..."
    "$SOD_BIN/sqssod" "${nxx}"
  fi
  echo " > Running gqssod for ${nxx}..."
  "$SOD_BIN/gqssod" "$nxx"
}

if ls -d n[0-9]*/ 2>/dev/null | grep -q .; then
  # Called from MAIN/: process all nXX/ folders
  check_main
  for nxx in $(ls -d n[0-9]*/ 2>/dev/null | sort); do
    run_for "${nxx%/}"
  done
elif [ -f ENSEMBLE ]; then
  # Called from inside nXX/: process this composition
  if [ ! -f ENERGIES ]; then
    echo "Error: ENERGIES not found in $(pwd)."
    exit 1
  fi
  # Check for required files in parent directory
  for f in EQMATRIX supercell.cif INSOD; do
    if [ ! -f "../$f" ]; then
      echo "Error: $f not found in parent directory."
      exit 1
    fi
  done
  # Check for TEMPERATURES in nXX/ or MAIN/
  if [ ! -f TEMPERATURES ] && [ ! -f ../TEMPERATURES ]; then
    echo "Error: TEMPERATURES not found in $(pwd) or in parent directory."
    exit 1
  fi
  nxx_name=$(basename "$(pwd)")
  if [ ! -f OUTSQS ]; then
    echo " > OUTSQS not found. Running sqssod first..."
    if [ ! -f INSQS ] && [ ! -f ../INSQS ]; then
      echo "Error: INSQS not found in $(pwd) or parent directory (MAIN/)."
      exit 1
    fi
    cd ..
    "$SOD_BIN/sqssod" "${nxx_name}"
    cd "${nxx_name}"
  fi
  # cd to parent (MAIN), then run gqssod with nXX arg
  cd ..
  echo " > Running gqssod for ${nxx_name}..."
  "$SOD_BIN/gqssod" "${nxx_name}"
else
  echo "Error: run sod_gqs.sh from MAIN/ or from an nXX/ folder."
  exit 1
fi
