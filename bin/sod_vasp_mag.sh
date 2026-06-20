#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${SCRIPT_DIR}/sod_common.sh"

SODPROJECT="$(sod_require_project_root "$PWD")" || exit 1
LEVEL_NAME="$(sod_find_enclosing_level_name "$SODPROJECT" "$PWD" || true)"

if [ -z "$LEVEL_NAME" ]; then
  # Called from SODPROJECT/: loop over all nXX/cYY/OUTCAR
  cd "$SODPROJECT" || exit 1
  if ! ls -d n[0-9]*/ 2>/dev/null | grep -q .; then
    echo "Error: no nXX/ folders found in SODPROJECT/."
    exit 1
  fi
  rm -f MAGNET
  for f in $(ls n*/c*/OUTCAR 2>/dev/null | sort); do
    grep "number of electron" "$f" | tail -1 | awk '{print $6}' >> MAGNET
  done
else
  # Called from nXX/: loop over cYY/OUTCAR in current folder
  cd "$SODPROJECT/$LEVEL_NAME" || exit 1
  if ! ls -d c[0-9]*/ 2>/dev/null | grep -q .; then
    echo "Error: no cYY/ folders found in ${LEVEL_NAME}/."
    exit 1
  fi
  rm -f MAGNET
  for f in $(ls c*/OUTCAR 2>/dev/null | sort); do
    grep "number of electron" "$f" | tail -1 | awk '{print $6}' >> MAGNET
  done
fi
