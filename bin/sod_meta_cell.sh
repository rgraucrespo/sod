#!/bin/bash

n_columns_ls=`ls -l |tail -1 |awk '{ FS = "|" } ; { print NF}'`
ls -l *.mout |awk -v nc=$n_columns_ls '{print "./getcell.sh " $nc}' > rungetcell
chmod +x rungetcell
rm -f a.dat b.dat c.dat
./rungetcell
paste  a.dat b.dat c.dat > cell.dat
awk '{print $1*$5, $1*$5*$9}' cell.dat > av.dat

