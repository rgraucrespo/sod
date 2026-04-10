#! /bin/bash

awk '$1 == "!" && $2 == "total" {printf "%.8f\n", $5 * 13.6056980659}' $1 | tail -1 >> ENERGIES
