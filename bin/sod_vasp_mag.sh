#!/bin/bash

rm -f MAGNET 

n_columns_ls=`ls -l |tail -1 |awk '{ FS = "|" } ; { print NF}'`
ls -l c*/OUTCAR |awk -v nc=$n_columns_ls '{print "magvasp.sh ", $nc}' > rungetmag
chmod +x rungetmag
./rungetmag
rm rungetmag

