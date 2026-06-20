#!/bin/bash
# This script generates the calculation input files (eg for VASP) from INSOD and ENSEMBLE

# Resolve script directory and find executables
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="${SCRIPT_DIR}:${PATH}"
. "${SCRIPT_DIR}/sod_common.sh"

SODPROJECT="$(sod_require_project_root "$PWD")" || exit 1
LAUNCH_DIR="$PWD"
LEVEL_NAME="$(sod_find_enclosing_level_name "$SODPROJECT" "$LAUNCH_DIR" || true)"
cd "$SODPROJECT" || exit 1

if [ -n "$LEVEL_NAME" ] && [ -f "$LAUNCH_DIR/ENSEMBLE" ]; then
  ENSEMBLE_SUBDIR="${LAUNCH_DIR#$SODPROJECT/}"
  SOD_ENSEMBLE_DIR="$ENSEMBLE_SUBDIR" genersod "$@"
else
  genersod "$@"
fi

# Read FILER value from last line of INSOD
FILER=$(tail -1 INSOD)

if [ $FILER -ne -1 ]; then

  if [ $FILER -eq 1 ]; then
    # GULP: genersod created n*/c*/ directories directly; just write a job_sender
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  (cd "$dir" && ${SOD_GULP:-gulp} < input.gin > output.gout)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ $FILER -eq 2 ]; then
    # LAMMPS: genersod created n*/c*/ directories directly; just write a job_sender
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  (cd "$dir" && ${SOD_LAMMPS:-lmp} < in.lammps > out.lammps)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ $FILER -eq 12 ]; then
    # CASTEP: genersod created n*/c*/ directories directly; just write a job_sender
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  (cd "$dir" && ${SOD_CASTEP:-castep} castep)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ $FILER -eq 11 ]; then
    # VASP: genersod created n*/c*/ directories directly; just write a job_sender
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  (cd "$dir" && ${SOD_VASP:-vasp})\ndone\n' > job_sender
    chmod +x job_sender

  elif [ $FILER -eq 13 ]; then
    # QE: genersod created n*/c*/ directories directly; just write a job_sender
    printf '#!/bin/bash\nfor dir in $(ls -d n*/c*/ 2>/dev/null | sort); do\n  (cd "$dir" && ${SOD_QE:-pw.x} < pw.in > pw.out)\ndone\n' > job_sender
    chmod +x job_sender

  elif [ $FILER -eq 0 ]; then
    # CIF: genersod creates n*/c*/configuration.cif directly; no calculation to run
    true

  else
    echo "Error: unknown FILER value $FILER"
    exit 1
  fi

fi
