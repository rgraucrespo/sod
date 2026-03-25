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


n_columns_ls=`ls -l |tail -1 |awk '{ FS = "|" } ; { print NF}'`
ls -l *.gout |awk -v nc=$n_columns_ls '{print "cellgulp.sh", $nc}' > rungetcell
chmod +x rungetcell
./rungetcell
paste  a.dat b.dat c.dat alpha.dat beta.dat gamma.dat volume.dat > cell.dat
rm rungetcell a.dat b.dat c.dat alpha.dat beta.dat gamma.dat volume.dat

# Read the first three letter of the argument and perform the calculations
argument=$(echo "$1" | cut -c1-3)
case "$argument" in
    cub)
        awk '{print $7^(1/3), $7 }' cell.dat > CELL
        ;;
    tet)
        awk '{print sqrt($7/$3), $3, $7 }' cell.dat > CELL
        ;;
    ort)
        awk '{print $7/($2*$3), $2, $3, $7 }' cell.dat > CELL
        ;;
    hex)
        awk '{print sqrt($7/(0.866*$3)), $3, $7 }' cell.dat > CELL
        ;;
    rho)
        awk '{print ($7/(sqrt(1-3*cos($4*3.141592/180)^2+2*cos($4*3.141592/180)^3)))^(1/3), $4, $7 }' cell.dat > CELL
        ;;
    mon)
        awk '{print $7/($2*$3*sin($5*3.141592/180)), $2, $3, $5, $7 }' cell.dat > CELL
        ;;
    tri)
        awk '{print $7/($2*$3*sqrt(1-cos($4*3.141592/180)^2-cos($5*3.141592/180)^2-cos($6*3.141592/180)^2+2*cos($4*3.141592/180)*cos($5*3.141592/180)*cos($6*3.141592/180))), $2, $3, $4, $5, $6, $7 }' cell.dat > CELL
        ;;
    *)
        echo "Error: non valid argument, it has to be one of the following: cub, tet, ort, hex, rho, mon, tri"
        exit 1
        ;;
esac


