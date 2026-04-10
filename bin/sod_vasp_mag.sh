#!/bin/bash

n_columns_ls=$(ls -l | tail -1 | awk '{ FS = "|" } ; { print NF}')

if ls -d n[0-9]*/ 2>/dev/null | grep -q .; then
  # Called from MAIN/: loop over all nXX/cYY/OUTCAR
  rm -f MAGNET
  ls -l n*/c*/OUTCAR | awk -v nc=$n_columns_ls '{print "magvasp.sh ", $nc}' > rungetmag
elif ls -d c[0-9]*/ 2>/dev/null | grep -q .; then
  # Called from nXX/: loop over cYY/OUTCAR in current folder
  rm -f MAGNET
  ls -l c*/OUTCAR | awk -v nc=$n_columns_ls '{print "magvasp.sh ", $nc}' > rungetmag
else
  echo "Error: run sod_vasp_mag.sh from MAIN/ or from an nXX/ folder."
  exit 1
fi

chmod +x rungetmag
./rungetmag
rm rungetmag
