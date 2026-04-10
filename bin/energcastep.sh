#! /bin/bash

awk '$1 == "Final" && $2 == "energy," {print $(NF-1)}' $1 | tail -1 >> ENERGIES
