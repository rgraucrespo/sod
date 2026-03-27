#! /bin/bash

awk     '{if (($1 == "Final") && ($2 == "free") && ($3 == "energy")) {print $5            }}' $1 >> ENERGIES 
