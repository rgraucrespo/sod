#!/bin/bash
# Run from SODPROJECT/ to process all nXX/ folders, or from any folder that
# contains ENSEMBLE and ENERGIES (e.g. nXX/ or nXX/PMEx/).

if ls -d n[0-9]*/ 2>/dev/null | grep -q .; then
  # Called from SODPROJECT/: run statsod in each nXX/ folder
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
elif [ -f ENSEMBLE ]; then
  # Called from a directory with ENSEMBLE (nXX/, nXX/PMEx/, etc.)
  statsod || { echo "Error: statsod failed"; exit 1; }
else
  echo "Error: run sod_stat.sh from SODPROJECT/ (to process all nXX/) or from any"
  echo "       directory containing ENSEMBLE and ENERGIES (e.g. nXX/ or nXX/PMEx/)."
  exit 1
fi
