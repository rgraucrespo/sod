#! /bin/bash

awk '{if (($1 == "Total") && ($3 == "energy") && ($6 == "eV") ) {print $5}}' $1 >> ENERGIES 
