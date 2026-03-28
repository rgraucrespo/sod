#!/bin/bash
# This script generates the calculation input files (eg for VASP) from INSOD and OUTSOD

genersod

FILER=$(<filer)

if [ $FILER -ne -1 ]; then

  if [ $FILER -eq 1 ]; then
    # GULP: genersod created n*/c*/ directories directly; just write a job_sender
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  (cd "$dir" && gulp < input.gin > output.gout)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ $FILER -eq 2 ]; then
    # LAMMPS: genersod created n*/c*/ directories directly; just write a job_sender
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  (cd "$dir" && lmp < in.lammps > out.lammps)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ $FILER -eq 12 ]; then
    # CASTEP: genersod created n*/c*/ directories directly; just write a job_sender
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  (cd "$dir" && castep castep)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ $FILER -eq 11 ]; then
    # VASP: genersod created n*/c*/ directories directly; just write a job_sender
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  (cd "$dir" && vasp_std)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ $FILER -eq 13 ]; then
    # QE: genersod created n*/c*/ directories directly; just write a job_sender
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  (cd "$dir" && pw.x < pw.in > pw.out)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ $FILER -eq 0 ]; then
    # CIF: genersod creates n*/c*/configuration.cif directly; no calculation to run
    true

  else
    echo "Error: unknown FILER value $FILER"
    exit 1
  fi

fi

rm -f filer
