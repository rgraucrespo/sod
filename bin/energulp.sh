#! /bin/bash

awk '{if (($1 == "Final") && ($2 == "energy") ) {print $4}}' $1 >> ENERGIES 
