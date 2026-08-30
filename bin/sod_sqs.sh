#!/bin/bash
# Run from SODPROJECT/, from an nXX/ folder, or from any subfolder containing ENSEMBLE and INSQS
# (e.g. nXX/random/).  EQMATRIX, supercell.cif, and INSOD must always be in SODPROJECT/.

SOD_BIN=$(cd "$(dirname "$0")" && pwd)
. "${SOD_BIN}/sod_common.sh"

SODPROJECT="$(sod_require_project_root "$PWD")" || exit 1
LAUNCH_DIR="$(cd "$PWD" && pwd -P)"
LEVEL_NAME="$(sod_find_enclosing_level_name "$SODPROJECT" "$LAUNCH_DIR" || true)"
cd "$SODPROJECT" || exit 1

check_sodproject() {
  if [ ! -f EQMATRIX ]; then
    echo "Error: EQMATRIX not found in $(pwd). Run combsod first."
    exit 1
  fi
}

# A sweep that scores nothing is an error, not a quiet success: the usual cause
# is a sampled ensemble sitting in nXX/random/ rather than nXX/, and the run
# would otherwise print a "Skipping" line and exit 0 having done nothing.
#
# Two counters, because the two ways of scoring nothing need different advice:
# `seen` counts directories that had an ENSEMBLE at all, `processed` counts the
# ones sqssod actually scored.  Reporting "no ENSEMBLE" when one was there and
# sqssod simply failed would send the user looking for the wrong problem.
seen=0
processed=0

run_for() {
  local ensemble_subdir="$1"
  if [ ! -f "${ensemble_subdir}/ENSEMBLE" ]; then
    echo " > Skipping ${ensemble_subdir}: ENSEMBLE not found."
    return 1
  fi
  seen=$((seen + 1))
  echo " > Running sqssod for ${ensemble_subdir}..."
  "$SOD_BIN/sqssod" "$ensemble_subdir" || return 1
  processed=$((processed + 1))
  return 0
}

# Nothing was scored. Say which of the two reasons it was, and fail.
nothing_scored_error() {
  if [ "$seen" -gt 0 ]; then
    echo "Error: sqssod failed for every ensemble it found; see its output above."
    exit 1
  fi
  echo "Error: no ENSEMBLE was found to score."
  local found
  found=$(find n[0-9]*/ -mindepth 2 -maxdepth 2 -type f -name ENSEMBLE 2>/dev/null | sort)
  if [ -n "$found" ]; then
    echo "       Ensembles were found in these subdirectories:"
    printf '%s\n' "$found" | sed 's|/ENSEMBLE$||; s|^|         |'
    echo "       sod_sqs.sh scores one level at a time; run it from one of them, e.g."
    echo "         cd $(printf '%s\n' "$found" | head -1 | sed 's|/ENSEMBLE$||') && sod_sqs.sh"
  else
    echo "       Run sod_comb.sh (enumeration) or sod_random.sh (sampling) first."
  fi
  exit 1
}

if [ -z "$LEVEL_NAME" ] && ls -d n[0-9]*/ 2>/dev/null | grep -q .; then
  # Called from SODPROJECT/: process all nXX/ folders
  check_sodproject
  for nxx in $(ls -d n[0-9]*/ 2>/dev/null | sort); do
    run_for "${nxx%/}" || true
  done
  if [ "$processed" -eq 0 ]; then
    nothing_scored_error
  fi
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
  run_for "$ENSEMBLE_SUBDIR" || true
  if [ "$processed" -eq 0 ]; then
    nothing_scored_error
  fi
else
  echo "Error: run sod_sqs.sh from SODPROJECT/ or from any SODPROJECT/nXX/ subfolder."
  exit 1
fi
