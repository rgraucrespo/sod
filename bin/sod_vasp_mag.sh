#!/bin/bash

if ls -d n[0-9]*/ 2>/dev/null | grep -q .; then
  # Called from MAIN/: loop over all nXX/cYY/OUTCAR
  rm -f MAGNET
  for f in $(ls n*/c*/OUTCAR 2>/dev/null | sort); do
    grep "number of electron" "$f" | tail -1 | awk '{print $6}' >> MAGNET
  done
elif ls -d c[0-9]*/ 2>/dev/null | grep -q .; then
  # Called from nXX/: loop over cYY/OUTCAR in current folder
  rm -f MAGNET
  for f in $(ls c*/OUTCAR 2>/dev/null | sort); do
    grep "number of electron" "$f" | tail -1 | awk '{print $6}' >> MAGNET
  done
else
  echo "Error: run sod_vasp_mag.sh from MAIN/ or from an nXX/ folder."
  exit 1
fi
