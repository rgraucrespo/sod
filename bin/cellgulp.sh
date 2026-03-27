#!/bin/bash

awk     '{if (($1 == "Final") && ($2 == "cell") ) {getline;getline;getline;print $2            }}' $1 >> a.dat
awk     '{if (($1 == "Final") && ($2 == "cell") ) {getline;getline;getline;getline; print $2         }}' $1 >> b.dat
awk     '{if (($1 == "Final") && ($2 == "cell") ) {getline;getline;getline;getline;getline; print $2 }}' $1 >> c.dat
awk     '{if (($1 == "Final") && ($2 == "cell") ) {getline;getline;getline;getline;getline; getline; print $2 }}' $1 >> alpha.dat
awk     '{if (($1 == "Final") && ($2 == "cell") ) {getline;getline;getline;getline;getline; getline; getline; print $2 }}' $1 >> beta.dat
awk     '{if (($1 == "Final") && ($2 == "cell") ) {getline;getline;getline;getline;getline; getline; getline; getline; print $2 }}' $1 >> gamma.dat
awk     '{if ($1 == "Non-primitive")  {print $5     }}' $1 >> volume.dat
