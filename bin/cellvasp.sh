#! /bin/bash

echo "final parameters" > tmp
head -5 $1 |tail -3 >> tmp

awk     '{if (($1 == "final") && ($2 == "parameters") ) {getline; print $0                  }}'  tmp >> a.dat
awk     '{if (($1 == "final") && ($2 == "parameters") ) {getline; getline; print $0         }}'  tmp >> b.dat
awk     '{if (($1 == "final") && ($2 == "parameters") ) {getline; getline; getline; print $0 }}' tmp >> c.dat

rm tmp

