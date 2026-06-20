#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${SCRIPT_DIR}/sod_common.sh"

SODPROJECT="$(sod_require_project_root "$PWD")" || exit 1
LEVEL_NAME="$(sod_find_enclosing_level_name "$SODPROJECT" "$PWD" || true)"

extract_energy() {
  local cdir="$1"
  awk '$1 == "!" && $2 == "total" {printf "%.8f\n", $5 * 13.6056980659}' "${cdir}pw.out" | tail -1
}

process_level() {
  local energies_file="$1"
  shift
  local cdirs=("$@")
  local n_missing=0
  rm -f "$energies_file"
  for cdir in "${cdirs[@]}"; do
    local m_raw
    m_raw=$(basename "${cdir%/}")
    local m=$((10#${m_raw#c}))
    if [ -f "${cdir}pw.out" ]; then
      local energy
      energy=$(extract_energy "$cdir")
      if [ -n "$energy" ]; then
        printf "%d  %s\n" "$m" "$energy" >> "$energies_file"
      else
        n_missing=$((n_missing + 1))
      fi
    else
      n_missing=$((n_missing + 1))
    fi
  done
  if [ "$n_missing" -gt 0 ]; then
    echo "Warning: missing energies for $n_missing configuration(s)."
  fi
}

if [ -z "$LEVEL_NAME" ]; then
  # Called from SODPROJECT/: extract energies for all nXX/ levels
  cd "$SODPROJECT" || exit 1
  if ! ls -d n[0-9]*/ 2>/dev/null | grep -q .; then
    echo "Error: no nXX/ folders found in SODPROJECT/."
    exit 1
  fi
  for ndir in $(ls -d n*/ 2>/dev/null | sort); do
    cdirs=()
    while IFS= read -r line; do
      cdirs+=("$line")
    done < <(ls -d "${ndir}"c*/ 2>/dev/null | sort)
    process_level "${ndir}ENERGIES" "${cdirs[@]}"
  done
else
  # Called from nXX/: extract energies for this level only
  cd "$SODPROJECT/$LEVEL_NAME" || exit 1
  if ! ls -d c[0-9]*/ 2>/dev/null | grep -q .; then
    echo "Error: no cYY/ folders found in ${LEVEL_NAME}/."
    exit 1
  fi
  cdirs=()
  while IFS= read -r line; do
    cdirs+=("$line")
  done < <(ls -d c*/ 2>/dev/null | sort)
  process_level "ENERGIES" "${cdirs[@]}"
fi
