#!/bin/bash
# This script generates the correct names for the input files
# It also generates a script that can be used to run all the simulations

# Resolve script directory and find executables
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="${SCRIPT_DIR}:${PATH}"
. "${SCRIPT_DIR}/sod_common.sh"

SODPROJECT="$(sod_require_project_root "$PWD")" || exit 1
cd "$SODPROJECT" || exit 1

rm -f ENSEMBLE supercell.cif EQMATRIX OPERATORS cSGO

clear

combsod || { echo "Error: combsod failed"; exit 1; }

# Read FILER value from last line of INSOD
if [ ! -f INSOD ]; then
  echo "Error: INSOD file not found"
  exit 1
fi
FILER=$(tail -1 INSOD)

if [ "$FILER" -ne -1 ]; then

  if [ "$FILER" -eq 1 ]; then
    # GULP: genersod creates n*/c*/ directories directly; just write a job_sender
    genersod || { echo "Error: genersod failed"; exit 1; }
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  echo " > Calculating: $dir"\n  (cd "$dir" && ${SOD_GULP:-gulp} < input.gin > output.gout)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ "$FILER" -eq 2 ]; then
    # LAMMPS: genersod creates n*/c*/ directories directly; just write a job_sender
    genersod || { echo "Error: genersod failed"; exit 1; }
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  echo " > Calculating: $dir"\n  (cd "$dir" && ${SOD_LAMMPS:-lmp} < in.lammps > out.lammps)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ "$FILER" -eq 12 ]; then
    # CASTEP: genersod creates n*/c*/ directories directly; just write a job_sender
    genersod || { echo "Error: genersod failed"; exit 1; }
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  echo " > Calculating: $dir"\n  (cd "$dir" && ${SOD_CASTEP:-castep} castep)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ "$FILER" -eq 11 ]; then
    # VASP: genersod creates n*/c*/ directories directly; just write a job_sender
    genersod || { echo "Error: genersod failed"; exit 1; }
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  echo " > Calculating: $dir"\n  (cd "$dir" && ${SOD_VASP:-vasp})\ndone\n' > job_sender
    chmod +x job_sender

  elif [ "$FILER" -eq 13 ]; then
    # QE: genersod creates n*/c*/ directories directly; just write a job_sender
    genersod || { echo "Error: genersod failed"; exit 1; }
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  echo " > Calculating: $dir"\n  (cd "$dir" && ${SOD_QE:-pw.x} < pw.in > pw.out)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ "$FILER" -eq 0 ]; then
    # CIF: genersod creates n*/c*/configuration.cif directly; no calculation to run
    genersod || { echo "Error: genersod failed"; exit 1; }

  else
    echo "Error: unknown FILER value $FILER"
    exit 1
  fi

fi
