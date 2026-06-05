#!/bin/bash
# Run from MAIN/ to score SQS for all nXX/ compositions, or from a single nXX/ folder.
# INSQS, EQMATRIX, supercell.cif, and INSOD must always be in MAIN/.

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
  echo " > Running sqssod for ${nxx}..."
  "$SOD_BIN/sqssod" "$nxx"
}

if ls -d n[0-9]*/ 2>/dev/null | grep -q .; then
  # Called from MAIN/: process all nXX/ folders
  check_main
  for nxx in $(ls -d n[0-9]*/ 2>/dev/null | sort); do
    run_for "${nxx%/}"
  done
elif [ -f ENSEMBLE ]; then
  # Called from inside nXX/: EQMATRIX etc. in MAIN/; INSQS in nXX/ (priority) or MAIN/
  for f in EQMATRIX supercell.cif INSOD; do
    if [ ! -f "../$f" ]; then
      echo "Error: $f not found in parent directory (MAIN/)."
      exit 1
    fi
  done
  if [ ! -f INSQS ] && [ ! -f ../INSQS ]; then
    echo "Error: INSQS not found in $(pwd) or parent directory (MAIN/)."
    exit 1
  fi
  nxx_name=$(basename "$(pwd)")
  cd ..
  echo " > Running sqssod for ${nxx_name}..."
  "$SOD_BIN/sqssod" "${nxx_name}"
else
  echo "Error: run sod_sqs.sh from MAIN/ or from an nXX/ folder."
  exit 1
fi
