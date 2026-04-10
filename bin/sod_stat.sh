#!/bin/bash
# Run from MAIN/ to process all nXX/ folders, or from a single nXX/ folder.

if ls -d n[0-9]*/ 2>/dev/null | grep -q .; then
  # Called from MAIN/: run statsod in each nXX/ folder
  first=1
  for ndir in $(ls -d n*/ 2>/dev/null | sort); do
    echo " > Running statsod in ${ndir}..."
    if [ $first -eq 1 ]; then
      (cd "$ndir" && statsod) || { echo "Error: statsod failed in ${ndir}"; exit 1; }
      first=0
    else
      (cd "$ndir" && statsod -q) || { echo "Error: statsod failed in ${ndir}"; exit 1; }
    fi
  done
elif [ -f OUTSOD ]; then
  # Called from nXX/: run statsod in the current folder
  statsod || { echo "Error: statsod failed"; exit 1; }
else
  echo "Error: run sod_stat.sh from MAIN/ or from an nXX/ folder."
  exit 1
fi
