#!/bin/bash
# Run from SODPROJECT/ to score GQS at all temperatures for all nXX/ compositions, or from a single nXX/ folder.
# Reads ENERGIES from nXX/ folder, TEMPERATURES from nXX/ (priority) or SODPROJECT/.
# Automatically runs sqssod first if OUTSQS is missing.

SOD_BIN=$(cd "$(dirname "$0")" && pwd)
. "${SOD_BIN}/sod_common.sh"

SODPROJECT="$(sod_require_project_root "$PWD")" || exit 1
LEVEL_NAME="$(sod_find_enclosing_level_name "$SODPROJECT" "$PWD" || true)"
cd "$SODPROJECT" || exit 1

check_sodproject() {
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

if [ -z "$LEVEL_NAME" ] && ls -d n[0-9]*/ 2>/dev/null | grep -q .; then
  # Called from SODPROJECT/: process all nXX/ folders
  check_sodproject
  for nxx in $(ls -d n[0-9]*/ 2>/dev/null | sort); do
    run_for "${nxx%/}"
  done
elif [ -n "$LEVEL_NAME" ]; then
  # Called from inside nXX/ or one of its subdirectories.
  if [ ! -f "${LEVEL_NAME}/ENERGIES" ]; then
    echo "Error: ENERGIES not found in ${LEVEL_NAME}/."
    exit 1
  fi
  for f in EQMATRIX supercell.cif INSOD; do
    if [ ! -f "$f" ]; then
      echo "Error: $f not found in SODPROJECT/."
      exit 1
    fi
  done
  # Check for TEMPERATURES in nXX/ or SODPROJECT/
  if [ ! -f "${LEVEL_NAME}/TEMPERATURES" ] && [ ! -f TEMPERATURES ]; then
    echo "Error: TEMPERATURES not found in ${LEVEL_NAME}/ or SODPROJECT/."
    exit 1
  fi
  if [ ! -f "${LEVEL_NAME}/OUTSQS" ]; then
    echo " > OUTSQS not found for ${LEVEL_NAME}. Running sqssod first..."
    if [ ! -f "${LEVEL_NAME}/INSQS" ] && [ ! -f INSQS ]; then
      echo "Error: INSQS not found in ${LEVEL_NAME}/ or SODPROJECT/."
      exit 1
    fi
    "$SOD_BIN/sqssod" "${LEVEL_NAME}"
  fi
  echo " > Running gqssod for ${LEVEL_NAME}..."
  "$SOD_BIN/gqssod" "${LEVEL_NAME}"
else
  echo "Error: run sod_gqs.sh from SODPROJECT/ or from any SODPROJECT/nXX/ subfolder."
  exit 1
fi
