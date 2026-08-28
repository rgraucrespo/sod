#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${SCRIPT_DIR}/sod_common.sh"

SODPROJECT="$(sod_require_project_root "$PWD")" || exit 1
LEVEL_NAME="$(sod_find_enclosing_level_name "$SODPROJECT" "$PWD" || true)"

extract_enthalpy() {
  local cdir="$1"
  awk '
    /energy\(sigma->0\)/ {
      for (i = 1; i <= NF; i++) {
        if ($i == "energy(sigma->0)") {
          j = i + 1
          if ($(j) == "=") j++
          energy = $(j)
        }
      }
    }
    /P[[:space:]]+V=/ {
      for (i = 1; i <= NF; i++) {
        if ($i == "V=") {
          pv = $(i + 1)
        }
      }
    }
    /P[[:space:]]+V=/ { have_pv = 1 }
    END {
      # No "P V=" line means no PSTRESS: an internal energy, not an enthalpy.
      if (energy != "" && have_pv) printf "%.10f\n", energy + pv
    }
  ' "${cdir}OUTCAR"
}

process_level() {
  local enthalpies_file="$1"
  local level_label="$2"
  shift 2
  local cdirs=("$@")
  local n_missing=0 n_outputs=0 n_written=0
  rm -f "$enthalpies_file"
  for cdir in "${cdirs[@]}"; do
    local m_raw
    m_raw=$(basename "${cdir%/}")
    local m=$((10#${m_raw#c}))
    if [ -f "${cdir}OUTCAR" ]; then
      n_outputs=$((n_outputs + 1))
      local enthalpy
      enthalpy=$(extract_enthalpy "$cdir")
      if [ -n "$enthalpy" ]; then
        printf "%d  %s\n" "$m" "$enthalpy" >> "$enthalpies_file"
        n_written=$((n_written + 1))
      else
        n_missing=$((n_missing + 1))
      fi
    else
      n_missing=$((n_missing + 1))
    fi
  done
  # With no "P V=" line anywhere the level was run at zero pressure, where
  # H = E exactly. Writing ENTHALPIES would only duplicate ENERGIES, so skip it
  # and clear any file left by an earlier constant-pressure run. This matches
  # sod_mace.sh, which writes ENTHALPIES only for nonzero -pressure.
  if [ "$n_written" -eq 0 ] && [ "$n_outputs" -gt 0 ]; then
    rm -f "$enthalpies_file"
    echo "Note: not writing ${level_label}ENTHALPIES: no "P V=" line found, so the level was run at zero pressure and H = E."
    return 0
  fi
  if [ "$n_missing" -gt 0 ]; then
    echo "Warning: missing enthalpies for $n_missing configuration(s)."
  fi
}

if [ -z "$LEVEL_NAME" ]; then
  # Called from SODPROJECT/: extract enthalpies for all nXX/ levels
  cd "$SODPROJECT" || exit 1
  if ! ls -d n[0-9]*/ 2>/dev/null | grep -q .; then
    echo "Error: no nXX/ folders found in SODPROJECT/."
    exit 1
  fi
  for ndir in $(ls -d n[0-9]*/ 2>/dev/null | sort); do
    cdirs=()
    while IFS= read -r name; do
      cdirs+=("${ndir}${name}/")
    done < <(sod_config_dirs_in_order "$ndir")
    process_level "${ndir}ENTHALPIES" "$ndir" "${cdirs[@]}"
  done
else
  # Called from nXX/: extract enthalpies for this level only
  cd "$SODPROJECT/$LEVEL_NAME" || exit 1
  if ! ls -d c[0-9]*/ 2>/dev/null | grep -q .; then
    echo "Error: no cYY/ folders found in ${LEVEL_NAME}/."
    exit 1
  fi
  cdirs=()
  while IFS= read -r name; do
    cdirs+=("${name}/")
  done < <(sod_config_dirs_in_order)
  process_level "ENTHALPIES" "" "${cdirs[@]}"
fi
