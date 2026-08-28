#!/bin/bash

# This script extracts unique cell parameters depending on the lattice system, and prints them in the CELL file.
# Usage: sod_vasp_cell.sh ARGUMENT
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
# CELL is positional, so it is written only when CONTCAR results cover every
# configuration listed by the level's ENSEMBLE.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${SCRIPT_DIR}/sod_common.sh"

SODPROJECT="$(sod_require_project_root "$PWD")" || exit 1
LEVEL_NAME="$(sod_find_enclosing_level_name "$SODPROJECT" "$PWD" || true)"

extract_vasp_cell() {
  local f="$1"
  printf 'final parameters\n' > _sod_cell_tmp
  head -5 "$f" | tail -3 >> _sod_cell_tmp
  awk '($1=="final") && ($2=="parameters") {getline; print $0}'             _sod_cell_tmp >> a.dat
  awk '($1=="final") && ($2=="parameters") {getline; getline; print $0}'    _sod_cell_tmp >> b.dat
  awk '($1=="final") && ($2=="parameters") {getline; getline; getline; print $0}' _sod_cell_tmp >> c.dat
  rm _sod_cell_tmp
}

# Read the first three letters of the argument and perform the calculations
argument=$(echo "$1" | cut -c1-3)
case "$argument" in
    cub|tet|ort|hex|rho|mon|tri) ;;
    *)
        echo "Error: non valid argument, it has to be one of the following: cub, tet, ort, hex, rho, mon, tri"
        exit 1
        ;;
esac

process_vasp_level() {
  local level_dir level_label rc cdir
  level_dir="$1"
  level_label="$(basename "$level_dir")"
  sod_level_outputs_complete "$level_dir" "CONTCAR"
  rc=$?
  if [ "$rc" -eq 2 ]; then
    # The ENSEMBLE is missing or unreadable, so coverage cannot be judged at
    # all. Leave any existing CELL alone and report a hard failure.
    echo "Error: not writing ${level_label}/CELL: its ENSEMBLE could not be read." >&2
    return 1
  fi
  if [ "$rc" -ne 0 ]; then
    rm -f "$level_dir/CELL"
    echo "Warning: not writing ${level_label}/CELL: CONTCAR results cover ${SOD_COMPLETE_OUTPUT_COUNT:-0} of ${SOD_ENSEMBLE_CONFIG_COUNT:-0} ENSEMBLE configurations, and CELL rows are positional." >&2
    return 0
  fi

  (
    cd "$level_dir" || exit 1
    rm -f cell.dat cellparams.dat a.dat b.dat c.dat
    while IFS= read -r cdir; do
      extract_vasp_cell "$cdir/CONTCAR"
    done < <(sod_config_dirs_in_order)
    paste a.dat b.dat c.dat > cell.dat
    rm a.dat b.dat c.dat

    # Convert raw Cartesian lattice vectors to a, b, c, alpha, beta, gamma, V.
    awk '{
      ax=$1; ay=$2; az=$3; bx=$4; by=$5; bz=$6; cx=$7; cy=$8; cz=$9
      a=sqrt(ax*ax+ay*ay+az*az); b=sqrt(bx*bx+by*by+bz*bz); c=sqrt(cx*cx+cy*cy+cz*cz)
      cosalpha=(bx*cx+by*cy+bz*cz)/(b*c); cosbeta=(ax*cx+ay*cy+az*cz)/(a*c)
      cosgamma=(ax*bx+ay*by+az*bz)/(a*b)
      alpha=atan2(sqrt(1-cosalpha^2),cosalpha)*180/3.14159265358979
      beta=atan2(sqrt(1-cosbeta^2),cosbeta)*180/3.14159265358979
      gamma=atan2(sqrt(1-cosgamma^2),cosgamma)*180/3.14159265358979
      V=ax*(by*cz-bz*cy)-ay*(bx*cz-bz*cx)+az*(bx*cy-by*cx); if(V<0)V=-V
      print a,b,c,alpha,beta,gamma,V
    }' cell.dat > cellparams.dat
    rm cell.dat

    case "$argument" in
      cub) awk '{print $7^(1/3),$7}' cellparams.dat > CELL.tmp ;;
      tet) awk '{print sqrt($7/$3),$3,$7}' cellparams.dat > CELL.tmp ;;
      ort) awk '{print $7/($2*$3),$2,$3,$7}' cellparams.dat > CELL.tmp ;;
      hex) awk '{print sqrt($7/(0.866*$3)),$3,$7}' cellparams.dat > CELL.tmp ;;
      rho) awk '{print ($7/sqrt(1-3*cos($4*3.141592/180)^2+2*cos($4*3.141592/180)^3))^(1/3),$4,$7}' cellparams.dat > CELL.tmp ;;
      mon) awk '{print $7/($2*$3*sin($5*3.141592/180)),$2,$3,$5,$7}' cellparams.dat > CELL.tmp ;;
      tri) awk '{print $7/($2*$3*sqrt(1-cos($4*3.141592/180)^2-cos($5*3.141592/180)^2-cos($6*3.141592/180)^2+2*cos($4*3.141592/180)*cos($5*3.141592/180)*cos($6*3.141592/180))),$2,$3,$4,$5,$6,$7}' cellparams.dat > CELL.tmp ;;
    esac
    rm cellparams.dat
    if [ "$(wc -l < CELL.tmp)" -ne "$SOD_ENSEMBLE_CONFIG_COUNT" ]; then
      rm -f CELL CELL.tmp
      echo "Error: not writing ${level_label}/CELL: one or more CONTCAR files contain no readable cell result." >&2
      exit 1
    fi
    mv CELL.tmp CELL
  ) || return 1
}

status=0
if [ -n "$LEVEL_NAME" ]; then
  process_vasp_level "$SODPROJECT/$LEVEL_NAME" || status=1
else
  found=0
  for level_dir in "$SODPROJECT"/n[0-9]*/; do
    [ -d "$level_dir" ] || continue
    found=1
    process_vasp_level "${level_dir%/}" || status=1
  done
  if [ "$found" -eq 0 ]; then
    echo "Error: no nXX/ folders found in SODPROJECT/." >&2
    exit 1
  fi
fi
exit "$status"
