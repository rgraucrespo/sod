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
    genersod || { echo "Error: genersod failed"; exit 1; }
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  (cd "$dir" && gulp < input.gin > output.gout)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ "$FILER" -eq 2 ]; then
    # LAMMPS: genersod creates n*/c*/ directories directly; just write a job_sender
    genersod || { echo "Error: genersod failed"; exit 1; }
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  (cd "$dir" && lmp < in.lammps > out.lammps)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ "$FILER" -eq 12 ]; then
    # CASTEP: genersod creates n*/c*/ directories directly; just write a job_sender
    genersod || { echo "Error: genersod failed"; exit 1; }
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  (cd "$dir" && castep castep)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ "$FILER" -eq 11 ]; then
    # VASP: genersod creates n*/c*/ directories directly; just write a job_sender
    genersod || { echo "Error: genersod failed"; exit 1; }
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  (cd "$dir" && vasp_std)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ "$FILER" -eq 13 ]; then
    # QE: genersod creates n*/c*/ directories directly; just write a job_sender
    genersod || { echo "Error: genersod failed"; exit 1; }
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  (cd "$dir" && pw.x < pw.in > pw.out)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ "$FILER" -eq 0 ]; then
    # CIF: genersod creates n*/c*/configuration.cif directly; no calculation to run
    genersod || { echo "Error: genersod failed"; exit 1; }

  else
    echo "Error: unknown FILER value $FILER"
    exit 1
  fi

fi

rm -f filer
