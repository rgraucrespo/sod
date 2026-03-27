#!/bin/bash
# This script generates the calculation input files (eg for VASP) from INSOD and OUTSOD

genersod

FILER=$(<filer)

if [ $FILER -ne 0 ]; then
  NSUBS=$(awk 'NR==1{print $1}' OUTSOD)
  FOLDER="n$(printf "%02d" $NSUBS)"
  if [ -d "$FOLDER" ]; then
    echo "Error: folder $FOLDER already exists. Remove it first to avoid overwriting."
    exit 1
  fi
  mkdir -p "$FOLDER"
  cd "$FOLDER"
  mv ../fort.* . 2>/dev/null || true
  cp ../OUTSOD .

  # FILER=1 for GULP
  if [ $FILER -eq 1 ];  then
    extin="gin"
    extout="gout"
    program="gulp"
  fi

  # FILER=2 for METADISE
  if [ $FILER -eq 2 ]; then
    extin="min"
    extout="mout"
    program="metadise"
  fi

  # FILER=11 for VASP
  if [ $FILER -eq 11 ]; then
    extin="vasp"
    extout="vout"
    program="vasp"
  fi

  # FILER=12 for CASTEP
  if [ $FILER -eq 12 ]; then
    extin="cell"
    extout="castep"
    program="castep"
  fi

  # FILER=13 for Quantum ESPRESSO
  if [ $FILER -eq 13 ]; then
    extin="pwi"
    extout="pwo"
    program="pw.x"
  fi

  ls fort.* > tmp1

  nconf=$(awk 'END{print NR}' tmp1)
  ndigits=${#nconf}
  sed s/fort.// tmp1 | awk -v nd="$ndigits" '{printf "%0" nd "d\n", $1-100000}' > tmp4

  awk -v extin=$extin '{print "c"$1"."extin}' tmp4 > tmp5

  paste tmp1 tmp5 > tmp6
  awk '{print "mv",$0}' tmp6 > tmp7
  chmod +x tmp7
  ./tmp7
  rm tmp*

  # This is to create the script that is going to run the calculations, with appropriate extension

 n_columns_ls=`ls -l |tail -1 |awk '{ FS = "|" } ; { print NF}'`
 ls -l *$extin |awk -v nc=$n_columns_ls -v extout=$extout -v program=$program '{print program,"<",$nc,">",$nc"."extout}' |sed s/$extin.$extout/$extout/  > job_sender
  chmod +x job_sender

  cd ..
fi

rm -f filer coord*
