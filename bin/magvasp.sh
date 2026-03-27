#! /bin/bash

grep "number of electron" $1|tail -1|awk '{print $6}' >> MAGNET
