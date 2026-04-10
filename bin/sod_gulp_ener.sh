#!/bin/bash

if ls -d n[0-9]*/ 2>/dev/null | grep -q .; then
  # Called from MAIN/: extract energies for all nXX/ levels
  for ndir in $(ls -d n*/ 2>/dev/null | sort); do
    rm -f "${ndir}ENERGIES"
    for cdir in $(ls -d "${ndir}"c*/ 2>/dev/null | sort); do
      if [ -f "${cdir}output.gout" ]; then
        awk '{if (($1 == "Final") && ($2 == "energy") ) {print $4}}' "${cdir}output.gout" >> "${ndir}ENERGIES"
      fi
    done
  done
elif ls -d c[0-9]*/ 2>/dev/null | grep -q .; then
  # Called from nXX/: extract energies for this level only
  rm -f ENERGIES
  for cdir in $(ls -d c*/ 2>/dev/null | sort); do
    if [ -f "${cdir}output.gout" ]; then
      awk '{if (($1 == "Final") && ($2 == "energy") ) {print $4}}' "${cdir}output.gout" >> ENERGIES
    fi
  done
else
  echo "Error: run sod_gulp_ener.sh from MAIN/ or from an nXX/ folder."
  exit 1
fi
