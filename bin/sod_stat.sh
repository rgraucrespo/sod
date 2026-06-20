#!/bin/bash
# Run from SODPROJECT/ to process all nXX/ folders, or from any folder that
# contains ENSEMBLE and ENERGIES (e.g. nXX/ or nXX/PMEx/).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="${SCRIPT_DIR}:${PATH}"
. "${SCRIPT_DIR}/sod_common.sh"

SODPROJECT="$(sod_find_project_root "$PWD" || true)"
LEVEL_NAME=""
if [ -n "$SODPROJECT" ]; then
  LEVEL_NAME="$(sod_find_enclosing_level_name "$SODPROJECT" "$PWD" || true)"
fi

if [ -f ENSEMBLE ]; then
  # Called from a directory with ENSEMBLE (nXX/, nXX/PMEx/, etc.)
  statsod || { echo "Error: statsod failed"; exit 1; }
elif [ -n "$SODPROJECT" ] && [ -n "$LEVEL_NAME" ] && [ -f "$SODPROJECT/$LEVEL_NAME/ENSEMBLE" ]; then
  # Called from below an nXX/ level, but not from the ENSEMBLE-bearing directory.
  (cd "$SODPROJECT/$LEVEL_NAME" && statsod) || { echo "Error: statsod failed in ${LEVEL_NAME}"; exit 1; }
elif [ -n "$SODPROJECT" ] && [ -z "$LEVEL_NAME" ] && ls -d "$SODPROJECT"/n[0-9]*/ 2>/dev/null | grep -q .; then
  # Called from SODPROJECT/: run statsod in each nXX/ folder
  cd "$SODPROJECT" || exit 1
  first=1
  for ndir in $(ls -d n*/ 2>/dev/null | sort); do
    echo " > Running statsod in ${ndir}..."
    if [ $first -eq 1 ]; then
      (cd "$ndir" && statsod) || { echo "Error: statsod failed in ${ndir}"; exit 1; }
      first=0
    else
      (cd "$ndir" && statsod -q) || { echo "Error: statsod failed in ${ndir}"; exit 1; }
    fi
  done
else
  echo "Error: run sod_stat.sh from SODPROJECT/ (to process all nXX/) or from any"
  echo "       directory containing ENSEMBLE and ENERGIES (e.g. nXX/ or nXX/PMEx/),"
  echo "       or from any subfolder below an nXX/ directory."
  exit 1
fi
