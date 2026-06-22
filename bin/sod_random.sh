#!/bin/bash
# sod_random.sh — uniform random sampling of the configuration space with SOD.
#
# Calls randomsod, which draws nconfigs independent uniform configurations at the
# INSOD target substitution level and writes them as nXX/random/ENSEMBLE. No
# energies are involved: the sample is Hamiltonian-independent, so energies (if
# wanted) are computed a posteriori via the normal structure-writer -> DFT ->
# statsod path, exactly as for an enumerated combsod ensemble. randomsod is the
# sampling counterpart of combsod, for substitution levels too large to enumerate.
#
# USAGE
#   sod_random.sh -nconf <N> [-sym on|off] [-seed clock|<int>]
#
#   -nconf <N>           Number of uniform draws (required). This is the number of
#                        draws, not the number of distinct configurations.
#   -sym on|off          Fold draws to symmetry representatives (default: off).
#                        With 'on', the degeneracy column holds visit counts, the
#                        correct importance weights for statsod canonical averages.
#   -seed clock|<int>    RNG seed (default: clock). Use a positive integer for a
#                        reproducible sample.
#
#   Run from the SODPROJECT directory, which must contain:
#     INSOD     — SOD input file (defines the target level and species)
#     SGO       — Space group operators file
#     EQMATRIX  — Equivalence matrix (only required with -symmetry on)
#
# OUTPUT
#     nXX/random/ENSEMBLE   (XX = target substitution level)

# Resolve script directory and find executables
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="${SCRIPT_DIR}:${PATH}"
. "${SCRIPT_DIR}/sod_common.sh"

SODPROJECT="$(sod_require_project_root "$PWD")" || exit 1
cd "$SODPROJECT" || exit 1

clear

# Parse arguments (defaults: sym on, seed clock). All flags are forwarded
# to randomsod, which re-validates them; we parse here for input-file checks.
SYMMETRY="off"
HAVE_NCONF=false
ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -nconf)
      [[ $# -ge 2 ]] || { echo "Error: -nconf requires a value."; exit 1; }
      ARGS+=("-nconf" "$2"); HAVE_NCONF=true; shift 2 ;;
    -sym)
      [[ $# -ge 2 ]] || { echo "Error: -sym requires a value (on|off)."; exit 1; }
      SYMMETRY="$2"; ARGS+=("-sym" "$2"); shift 2 ;;
    -seed)
      [[ $# -ge 2 ]] || { echo "Error: -seed requires a value (clock|<int>)."; exit 1; }
      ARGS+=("-seed" "$2"); shift 2 ;;
    *)
      echo "Error: unrecognised argument to sod_random.sh: $1"
      echo "Usage: sod_random.sh -nconf <N> [-sym on|off] [-seed clock|<int>]"
      exit 1 ;;
  esac
done

if ! $HAVE_NCONF; then
  echo "Error: -nconf is required."
  echo "Usage: sod_random.sh -nconf <N> [-sym on|off] [-seed clock|<int>]"
  exit 1
fi

# Required input files: INSOD and SGO always; EQMATRIX only when reducing.
for reqfile in INSOD SGO; do
  if [ ! -f "$reqfile" ]; then
    echo "Error: required file '$reqfile' not found in SODPROJECT/."
    exit 1
  fi
done

if [ "$SYMMETRY" = "on" ] && [ ! -f EQMATRIX ]; then
  echo "Error: EQMATRIX not found, required for -sym on."
  echo "Generate it with sod_comb.sh, or run with -sym off."
  exit 1
fi

randomsod "${ARGS[@]}" || exit 1

exit 0
