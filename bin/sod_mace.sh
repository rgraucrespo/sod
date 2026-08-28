#!/bin/bash
# sod_mace.sh — MACE machine-learning-potential energies and relaxation.
#
# Wrapper around pysod/sod_mace.py. Like the other post-processing scripts it
# runs from SODPROJECT/ (every nXX/) or from SODPROJECT/nXX/ (that level only),
# and it reads options from a mace_settings.yaml file when one is present
# (nXX/ takes priority over SODPROJECT). Anything given on the command line
# overrides the file.
#
# The Python stack (torch, ase, mace-torch, nvalchemi-toolkit) is optional and
# is not installed by 'make'. Set SOD_PYTHON to the interpreter that has it:
#
#   export SOD_PYTHON=~/miniconda3/envs/nvalchemi/bin/python
#
# See pysod/README.md and docs/mace.rst.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${SCRIPT_DIR}/sod_common.sh"

# Fail early with the standard message if we are not inside a SODPROJECT; the
# Python tool checks this too, but this keeps the error identical to the other
# wrappers and avoids paying interpreter start-up to learn it.
sod_require_project_root "$PWD" >/dev/null || exit 1

# Locate the Python implementation: source tree first, then an installed tree
# (make install puts bin/*.sh in PREFIX/bin and pysod in PREFIX/lib/sod/pysod).
SOD_MACE_PY=""
for candidate in \
    "${SCRIPT_DIR}/../pysod/sod_mace.py" \
    "${SCRIPT_DIR}/../lib/sod/pysod/sod_mace.py"
do
    if [ -f "$candidate" ]; then
        SOD_MACE_PY="$candidate"
        break
    fi
done
if [ -z "$SOD_MACE_PY" ]; then
    echo "Error: could not locate pysod/sod_mace.py relative to ${SCRIPT_DIR}." >&2
    echo "       Expected ../pysod/ (source tree) or ../lib/sod/pysod/ (installed)." >&2
    exit 1
fi

PY="${SOD_PYTHON:-python3}"
if ! command -v "$PY" >/dev/null 2>&1; then
    echo "Error: python interpreter not found: $PY" >&2
    echo "       Set SOD_PYTHON to an interpreter with the MACE stack installed," >&2
    echo "       for example: export SOD_PYTHON=~/miniconda3/envs/nvalchemi/bin/python" >&2
    exit 1
fi
if ! "$PY" -c 'import torch, ase, mace, nvalchemi' >/dev/null 2>&1; then
    echo "Error: the MACE stack is not available in $PY." >&2
    echo "       Needs torch, ase, mace-torch and nvalchemi-toolkit." >&2
    echo "       Set SOD_PYTHON to the right interpreter, or see pysod/README.md" >&2
    echo "       for installation instructions." >&2
    exit 1
fi

exec "$PY" "$SOD_MACE_PY" "$@"
