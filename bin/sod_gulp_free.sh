#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${SCRIPT_DIR}/sod_common.sh"

SODPROJECT="$(sod_require_project_root "$PWD")" || exit 1
LEVEL_NAME="$(sod_find_enclosing_level_name "$SODPROJECT" "$PWD" || true)"

if [ -z "$LEVEL_NAME" ]; then
  # Called from SODPROJECT/: extract free energies for all nXX/ levels
  cd "$SODPROJECT" || exit 1
  if ! ls -d n[0-9]*/ 2>/dev/null | grep -q .; then
    echo "Error: no nXX/ folders found in SODPROJECT/."
    exit 1
  fi
  for ndir in $(ls -d n*/ 2>/dev/null | sort); do
    rm -f "${ndir}ENERGIES"
    for cdir in $(ls -d "${ndir}"c*/ 2>/dev/null | sort); do
      if [ -f "${cdir}output.gout" ]; then
        awk '($1=="Final") && ($2=="free") && ($3=="energy") {print $5}' "${cdir}output.gout" >> ENERGIES
      fi
    done
  done
else
  # Called from nXX/: extract free energies for this level only
  cd "$SODPROJECT/$LEVEL_NAME" || exit 1
  if ! ls -d c[0-9]*/ 2>/dev/null | grep -q .; then
    echo "Error: no cYY/ folders found in ${LEVEL_NAME}/."
    exit 1
  fi
  rm -f ENERGIES
  for cdir in $(ls -d c*/ 2>/dev/null | sort); do
    if [ -f "${cdir}output.gout" ]; then
      awk '($1=="Final") && ($2=="free") && ($3=="energy") {print $5}' "${cdir}output.gout" >> ENERGIES
    fi
  done
fi
