#!/bin/bash
# This script generates the correct names for the input files
# It also generates a script that can be used to run all the simulations

rm -f OUTSOD supercell.cif EQMATRIX OPERATORS cSGO

clear

combsod || { echo "Error: combsod failed"; exit 1; }

if [ ! -f filer ]; then
  echo "Error: filer not created by combsod"
  exit 1
fi
FILER=$(<filer)

if [ "$FILER" -ne -1 ]; then

  if [ "$FILER" -eq 1 ]; then
    # GULP: genersod creates n*/c*/ directories directly; just write a job_sender
    genersod
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  (cd "$dir" && gulp < input.gin > output.gout)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ "$FILER" -eq 2 ]; then
    # LAMMPS: genersod creates n*/c*/ directories directly; just write a job_sender
    genersod
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  (cd "$dir" && lmp < in.lammps > out.lammps)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ "$FILER" -eq 12 ]; then
    # CASTEP: genersod creates n*/c*/ directories directly; just write a job_sender
    genersod
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  (cd "$dir" && castep castep)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ "$FILER" -eq 11 ]; then
    # VASP: genersod creates n*/c*/ directories directly; just write a job_sender
    genersod
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  (cd "$dir" && vasp_std)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ "$FILER" -eq 13 ]; then
    # QE: genersod creates n*/c*/ directories directly; just write a job_sender
    genersod
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  (cd "$dir" && pw.x < pw.in > pw.out)\ndone\n' > job_sender
    chmod +x job_sender

  else
    NSUBS=$(awk 'NR==1{print $1}' OUTSOD)
    FOLDER="n$(printf "%02d" $NSUBS)"
    if [ -d "$FOLDER" ]; then
      echo "Error: folder $FOLDER already exists. Remove it first to avoid overwriting."
      exit 1
    fi
    genersod
    mkdir -p "$FOLDER"
    cd "$FOLDER"
    mv ../fort.* . 2>/dev/null || true
    cp ../OUTSOD .

    # FILER=0 for CIF
    if [ $FILER -eq 0 ]; then
      extin="cif"
      extout="cout"
      program="none"
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

fi

rm -f filer
