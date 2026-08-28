#!/bin/bash

# This script extracts unique cell parameters depending on the lattice system, and prints them in the CELL file.
# Usage: sod_gulp_cell.sh ARGUMENT
# ARGUMENT must be one of the following cases: "cubic", "tetragonal", "orthorhombic", "hexagonal", "rhombohedral", "monoclinic" or "triclinic"
# (it is enough to specify the first three letters, e.g. "cub")
# The cell parameters written in each case are as follows:
#
#cub  a V
#tet  a c V
#ort  a b c V
#hex  a c V
#rho  a alpha V
#mon  a b c beta V
#tri  a b c alpha beta gamma V
#
# CELL is positional, so it is written only when output.gout results cover every
# configuration listed by the level's ENSEMBLE.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${SCRIPT_DIR}/sod_common.sh"

SODPROJECT="$(sod_require_project_root "$PWD")" || exit 1
LEVEL_NAME="$(sod_find_enclosing_level_name "$SODPROJECT" "$PWD" || true)"

# Read the first three letter of the argument and perform the calculations
argument=$(echo "$1" | cut -c1-3)
case "$argument" in
    cub|tet|ort|hex|rho|mon|tri) ;;
    *)
        echo "Error: non valid argument, it has to be one of the following: cub, tet, ort, hex, rho, mon, tri"
        exit 1
        ;;
esac

process_gulp_level() {
  local level_dir level_label rc cdir dir
  level_dir="$1"
  level_label="$(basename "$level_dir")"
  sod_level_outputs_complete "$level_dir" "output.gout"
  rc=$?
  if [ "$rc" -eq 2 ]; then
    # The ENSEMBLE is missing or unreadable, so coverage cannot be judged at
    # all. Leave any existing CELL alone and report a hard failure.
    echo "Error: not writing ${level_label}/CELL: its ENSEMBLE could not be read." >&2
    return 1
  fi
  if [ "$rc" -ne 0 ]; then
    rm -f "$level_dir/CELL"
    echo "Warning: not writing ${level_label}/CELL: GULP results cover ${SOD_COMPLETE_OUTPUT_COUNT:-0} of ${SOD_ENSEMBLE_CONFIG_COUNT:-0} ENSEMBLE configurations, and CELL rows are positional." >&2
    return 0
  fi

  (
    cd "$level_dir" || exit 1
    rm -f cell.dat a.dat b.dat c.dat alpha.dat beta.dat gamma.dat volume.dat
    while IFS= read -r cdir; do
      # Coverage is already complete, so a cYY without output.gout is a spare
      # directory rather than a gap. Skip it instead of letting awk report it.
      [ -f "$cdir/output.gout" ] || continue
      dir="$cdir/"
      awk '($1=="Final") && ($2=="cell") {getline;getline;getline;             print $2}' "${dir}output.gout" >> a.dat
      awk '($1=="Final") && ($2=="cell") {getline;getline;getline;getline;     print $2}' "${dir}output.gout" >> b.dat
      awk '($1=="Final") && ($2=="cell") {getline;getline;getline;getline;getline; print $2}' "${dir}output.gout" >> c.dat
      awk '($1=="Final") && ($2=="cell") {getline;getline;getline;getline;getline;getline; print $2}' "${dir}output.gout" >> alpha.dat
      awk '($1=="Final") && ($2=="cell") {getline;getline;getline;getline;getline;getline;getline; print $2}' "${dir}output.gout" >> beta.dat
      awk '($1=="Final") && ($2=="cell") {getline;getline;getline;getline;getline;getline;getline;getline; print $2}' "${dir}output.gout" >> gamma.dat
      awk '$1=="Non-primitive" {print $5}' "${dir}output.gout" >> volume.dat
    done < <(sod_config_dirs_in_order)
    paste a.dat b.dat c.dat alpha.dat beta.dat gamma.dat volume.dat > cell.dat
    rm a.dat b.dat c.dat alpha.dat beta.dat gamma.dat volume.dat

    case "$argument" in
      cub) awk '{print $7^(1/3),$7}' cell.dat > CELL.tmp ;;
      tet) awk '{print sqrt($7/$3),$3,$7}' cell.dat > CELL.tmp ;;
      ort) awk '{print $7/($2*$3),$2,$3,$7}' cell.dat > CELL.tmp ;;
      hex) awk '{print sqrt($7/(0.866*$3)),$3,$7}' cell.dat > CELL.tmp ;;
      rho) awk '{print ($7/sqrt(1-3*cos($4*3.141592/180)^2+2*cos($4*3.141592/180)^3))^(1/3),$4,$7}' cell.dat > CELL.tmp ;;
      mon) awk '{print $7/($2*$3*sin($5*3.141592/180)),$2,$3,$5,$7}' cell.dat > CELL.tmp ;;
      tri) awk '{print $7/($2*$3*sqrt(1-cos($4*3.141592/180)^2-cos($5*3.141592/180)^2-cos($6*3.141592/180)^2+2*cos($4*3.141592/180)*cos($5*3.141592/180)*cos($6*3.141592/180))),$2,$3,$4,$5,$6,$7}' cell.dat > CELL.tmp ;;
    esac
    rm cell.dat
    if [ "$(wc -l < CELL.tmp)" -ne "$SOD_ENSEMBLE_CONFIG_COUNT" ]; then
      rm -f CELL CELL.tmp
      echo "Error: not writing ${level_label}/CELL: one or more output.gout files contain no readable cell result." >&2
      exit 1
    fi
    mv CELL.tmp CELL
  ) || return 1
}

status=0
if [ -n "$LEVEL_NAME" ]; then
  process_gulp_level "$SODPROJECT/$LEVEL_NAME" || status=1
else
  found=0
  for level_dir in "$SODPROJECT"/n[0-9]*/; do
    [ -d "$level_dir" ] || continue
    found=1
    process_gulp_level "${level_dir%/}" || status=1
  done
  if [ "$found" -eq 0 ]; then
    echo "Error: no nXX/ folders found in SODPROJECT/." >&2
    exit 1
  fi
fi
exit "$status"
