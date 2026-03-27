#! /bin/bash

grep sigma $1 |awk '{print $7}' |tail -1 >> ENERGIES

