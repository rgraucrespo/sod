#!/bin/bash

rm -f ENERGIES 
n_columns_ls=`ls -l |tail -1 |awk '{ FS = "|" } ; { print NF}'`
ls -l */OUTCAR |awk -v nc=$n_columns_ls '{print "enervasp.sh ", $nc}'> rungetener
chmod +x rungetener
./rungetener
rm rungetener

