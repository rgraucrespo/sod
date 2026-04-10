#!/bin/bash

if ls -d n[0-9]*/ 2>/dev/null | grep -q .; then
  # Called from MAIN/: extract energies for all nXX/ levels
  for ndir in $(ls -d n*/ 2>/dev/null | sort); do
    rm -f "${ndir}ENERGIES"
    for cdir in $(ls -d "${ndir}"c*/ 2>/dev/null | sort); do
      if [ -f "${cdir}pw.out" ]; then
        awk '$1 == "!" && $2 == "total" {printf "%.8f\n", $5 * 13.6056980659}' "${cdir}pw.out" | tail -1 >> "${ndir}ENERGIES"
      fi
    done
  done
elif ls -d c[0-9]*/ 2>/dev/null | grep -q .; then
  # Called from nXX/: extract energies for this level only
  rm -f ENERGIES
  for cdir in $(ls -d c*/ 2>/dev/null | sort); do
    if [ -f "${cdir}pw.out" ]; then
      awk '$1 == "!" && $2 == "total" {printf "%.8f\n", $5 * 13.6056980659}' "${cdir}pw.out" | tail -1 >> ENERGIES
    fi
  done
else
  echo "Error: run sod_qe_ener.sh from MAIN/ or from an nXX/ folder."
  exit 1
fi
