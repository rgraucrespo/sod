#!/bin/bash

rm -f ENERGIES

for dir in $(ls -d n*/c*/ 2>/dev/null | sort); do
  if [ -f "${dir}output.gout" ]; then
    energulp.sh "${dir}output.gout"
  fi
done
