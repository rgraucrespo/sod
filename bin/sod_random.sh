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
#   sod_random.sh -nconfigs <N> [-symmetry on|off] [-seed clock|<int>]
#
#   -nconfigs <N>        Number of uniform draws (required). This is the number of
#                        draws, not the number of distinct configurations.
#   -symmetry on|off     Fold draws to symmetry representatives (default: on).
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

# Parse arguments (defaults: symmetry on, seed clock). All flags are forwarded
# to randomsod, which re-validates them; we parse here for input-file checks.
SYMMETRY="on"
ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -nconfigs)
      [[ $# -ge 2 ]] || { echo "Error: -nconfigs requires a value."; exit 1; }
      ARGS+=("-nconfigs" "$2"); shift 2 ;;
    -symmetry)
      [[ $# -ge 2 ]] || { echo "Error: -symmetry requires a value (on|off)."; exit 1; }
      SYMMETRY="$2"; ARGS+=("-symmetry" "$2"); shift 2 ;;
    -seed)
      [[ $# -ge 2 ]] || { echo "Error: -seed requires a value (clock|<int>)."; exit 1; }
      ARGS+=("-seed" "$2"); shift 2 ;;
    *)
      echo "Error: unrecognised argument to sod_random.sh: $1"
      echo "Usage: sod_random.sh -nconfigs <N> [-symmetry on|off] [-seed clock|<int>]"
      exit 1 ;;
  esac
done

# Required input files: INSOD and SGO always; EQMATRIX only when reducing.
for reqfile in INSOD SGO; do
  if [ ! -f "$reqfile" ]; then
    echo "Error: required file '$reqfile' not found in SODPROJECT/."
    exit 1
  fi
done

if [ "$SYMMETRY" = "on" ] && [ ! -f EQMATRIX ]; then
  echo "Error: EQMATRIX not found, required for -symmetry on."
  echo "Generate it with sod_comb.sh, or run with -symmetry off."
  exit 1
fi

randomsod "${ARGS[@]}" || exit 1

exit 0
