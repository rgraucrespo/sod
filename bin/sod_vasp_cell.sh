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

rm -f cell.dat a.dat b.dat c.dat

n_columns_ls=`ls -l |tail -1 |awk '{ FS = "|" } ; { print NF}'`
ls -l */CONTCAR |awk -v nc=$n_columns_ls '{print "cellvasp.sh", $nc}' > rungetcell
chmod +x rungetcell
./rungetcell
paste a.dat b.dat c.dat > cell.dat
rm rungetcell a.dat b.dat c.dat

# Convert raw Cartesian lattice vectors to a, b, c, alpha, beta, gamma, V
# cell.dat columns: ax ay az bx by bz cx cy cz
awk '{
  ax=$1; ay=$2; az=$3
  bx=$4; by=$5; bz=$6
  cx=$7; cy=$8; cz=$9
  a = sqrt(ax*ax + ay*ay + az*az)
  b = sqrt(bx*bx + by*by + bz*bz)
  c = sqrt(cx*cx + cy*cy + cz*cz)
  cosalpha = (bx*cx + by*cy + bz*cz) / (b*c)
  cosbeta  = (ax*cx + ay*cy + az*cz) / (a*c)
  cosgamma = (ax*bx + ay*by + az*bz) / (a*b)
  alpha = atan2(sqrt(1-cosalpha^2), cosalpha) * 180 / 3.14159265358979
  beta  = atan2(sqrt(1-cosbeta^2),  cosbeta)  * 180 / 3.14159265358979
  gamma = atan2(sqrt(1-cosgamma^2), cosgamma) * 180 / 3.14159265358979
  V = ax*(by*cz - bz*cy) - ay*(bx*cz - bz*cx) + az*(bx*cy - by*cx)
  if (V < 0) V = -V
  print a, b, c, alpha, beta, gamma, V
}' cell.dat > cellparams.dat
rm cell.dat

# Read the first three letters of the argument and perform the calculations
argument=$(echo "$1" | cut -c1-3)
case "$argument" in
    cub)
        awk '{print $7^(1/3), $7 }' cellparams.dat > CELL
        ;;
    tet)
        awk '{print sqrt($7/$3), $3, $7 }' cellparams.dat > CELL
        ;;
    ort)
        awk '{print $7/($2*$3), $2, $3, $7 }' cellparams.dat > CELL
        ;;
    hex)
        awk '{print sqrt($7/(0.866*$3)), $3, $7 }' cellparams.dat > CELL
        ;;
    rho)
        awk '{print ($7/(sqrt(1-3*cos($4*3.141592/180)^2+2*cos($4*3.141592/180)^3)))^(1/3), $4, $7 }' cellparams.dat > CELL
        ;;
    mon)
        awk '{print $7/($2*$3*sin($5*3.141592/180)), $2, $3, $5, $7 }' cellparams.dat > CELL
        ;;
    tri)
        awk '{print $7/($2*$3*sqrt(1-cos($4*3.141592/180)^2-cos($5*3.141592/180)^2-cos($6*3.141592/180)^2+2*cos($4*3.141592/180)*cos($5*3.141592/180)*cos($6*3.141592/180))), $2, $3, $4, $5, $6, $7 }' cellparams.dat > CELL
        ;;
    *)
        echo "Error: non valid argument, it has to be one of the following: cub, tet, ort, hex, rho, mon, tri"
        exit 1
        ;;
esac

rm cellparams.dat
