#!/bin/bash
# Run from SODPROJECT/, from an nXX/ folder, or from any subfolder containing ENSEMBLE and INSQS
# (e.g. nXX/random/).  EQMATRIX, supercell.cif, and INSOD must always be in SODPROJECT/.

SOD_BIN=$(cd "$(dirname "$0")" && pwd)
. "${SOD_BIN}/sod_common.sh"

SODPROJECT="$(sod_require_project_root "$PWD")" || exit 1
LAUNCH_DIR="$PWD"
LEVEL_NAME="$(sod_find_enclosing_level_name "$SODPROJECT" "$LAUNCH_DIR" || true)"
cd "$SODPROJECT" || exit 1

check_sodproject() {
  if [ ! -f EQMATRIX ]; then
    echo "Error: EQMATRIX not found in $(pwd). Run combsod first."
    exit 1
  fi
}

run_for() {
  local ensemble_subdir="$1"
  if [ ! -f "${ensemble_subdir}/ENSEMBLE" ]; then
    echo " > Skipping ${ensemble_subdir}: ENSEMBLE not found."
    return
  fi
  echo " > Running sqssod for ${ensemble_subdir}..."
  "$SOD_BIN/sqssod" "$ensemble_subdir"
}

if [ -z "$LEVEL_NAME" ] && ls -d n[0-9]*/ 2>/dev/null | grep -q .; then
  # Called from SODPROJECT/: process all nXX/ folders
  check_sodproject
  for nxx in $(ls -d n[0-9]*/ 2>/dev/null | sort); do
    run_for "${nxx%/}"
  done
elif [ -n "$LEVEL_NAME" ]; then
  # Called from inside nXX/ or one of its subdirectories (e.g. nXX/random/).
  # Pass the path relative to SODPROJECT/ so sqssod reads ENSEMBLE and INSQS from there.
  ENSEMBLE_SUBDIR="${LAUNCH_DIR#$SODPROJECT/}"
  for f in EQMATRIX supercell.cif INSOD; do
    if [ ! -f "$f" ]; then
      echo "Error: $f not found in SODPROJECT/."
      exit 1
    fi
  done
  if [ ! -f "${ENSEMBLE_SUBDIR}/INSQS" ] && [ ! -f INSQS ]; then
    echo "Error: INSQS not found in ${ENSEMBLE_SUBDIR}/ or SODPROJECT/."
    exit 1
  fi
  run_for "$ENSEMBLE_SUBDIR"
else
  echo "Error: run sod_sqs.sh from SODPROJECT/ or from any SODPROJECT/nXX/ subfolder."
  exit 1
fi
