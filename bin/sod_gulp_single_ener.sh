#!/bin/bash

extract_energy() {
  local cdir="$1"
  awk '($1=="Total") && ($3=="energy") && ($6=="eV") {print $5}' "${cdir}output.gout"
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
    if [ -f "${cdir}output.gout" ]; then
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

if ls -d n[0-9]*/ 2>/dev/null | grep -q .; then
  # Called from MAIN/: extract single-point energies for all nXX/ levels
  for ndir in $(ls -d n*/ 2>/dev/null | sort); do
    cdirs=()
    while IFS= read -r line; do
      cdirs+=("$line")
    done < <(ls -d "${ndir}"c*/ 2>/dev/null | sort)
    process_level "${ndir}ENERGIES" "${cdirs[@]}"
  done
elif ls -d c[0-9]*/ 2>/dev/null | grep -q .; then
  # Called from nXX/: extract single-point energies for this level only
  cdirs=()
  while IFS= read -r line; do
    cdirs+=("$line")
  done < <(ls -d c*/ 2>/dev/null | sort)
  process_level "ENERGIES" "${cdirs[@]}"
else
  echo "Error: run sod_gulp_single_ener.sh from MAIN/ or from an nXX/ folder."
  exit 1
fi
