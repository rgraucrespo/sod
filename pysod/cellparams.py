"""Lattice parameters and the SOD ``CELL`` file.

Python port of the geometry in ``bin/sod_vasp_cell.sh`` / ``bin/sod_gulp_cell.sh``.
Those scripts turn each relaxed configuration's lattice vectors into a row of
``nXX/CELL``, with the columns chosen by the lattice system, for averaging over
the ensemble with ``statsod``.

Two quirks of those scripts are reproduced deliberately, so that a CELL written
here is directly comparable with one written by the extractors:

1. The first parameter is back-derived from the volume rather than measured
   (cubic ``a`` is ``V**(1/3)``, not ``|a|``) -- it is the volume-equivalent
   lattice parameter.
2. The shell pipeline passes ``a b c alpha beta gamma V`` through awk's default
   ``OFMT`` (``%.6g``) into its intermediate ``cellparams.dat``, so the final
   columns are computed from **already-rounded** values.  Keeping full precision
   instead shifts the sixth significant figure of ``a`` for about one row in
   eight (~1e-5 A, physically irrelevant, but a visible diff).

Both are verified against ``examples/example01/FILER11_vasp/n04/CELL``, which
``sod_vasp_cell.sh cub`` produced from the committed DFT CONTCARs: this module
reproduces all 71 rows exactly.
"""

from __future__ import annotations

import math

# Column layout of CELL for each lattice system, as documented in
# bin/sod_vasp_cell.sh.
CELL_COLUMNS = {
    "cub": "a V",
    "tet": "a c V",
    "ort": "a b c V",
    "hex": "a c V",
    "rho": "a alpha V",
    "mon": "a b c beta V",
    "tri": "a b c alpha beta gamma V",
}

LATTICE_SYSTEMS = tuple(CELL_COLUMNS)


class LatticeSystemError(ValueError):
    """Raised for an unrecognised lattice system."""


def normalise_system(name: str) -> str:
    """Accept the first three letters, as the shell scripts do."""
    key = str(name).strip().lower()[:3]
    if key not in CELL_COLUMNS:
        raise LatticeSystemError(
            f"non valid lattice system: {name}. "
            f"It has to be one of: {', '.join(LATTICE_SYSTEMS)}"
        )
    return key


def cell_parameters(cell) -> tuple[float, float, float, float, float, float, float]:
    """Lattice vectors (3x3, vectors as rows) -> a, b, c, alpha, beta, gamma, V."""
    (ax, ay, az), (bx, by, bz), (cx, cy, cz) = (
        (float(cell[i][0]), float(cell[i][1]), float(cell[i][2])) for i in range(3)
    )
    a = math.sqrt(ax * ax + ay * ay + az * az)
    b = math.sqrt(bx * bx + by * by + bz * bz)
    c = math.sqrt(cx * cx + cy * cy + cz * cz)
    cosalpha = (bx * cx + by * cy + bz * cz) / (b * c)
    cosbeta = (ax * cx + ay * cy + az * cz) / (a * c)
    cosgamma = (ax * bx + ay * by + az * bz) / (a * b)
    alpha = math.degrees(math.atan2(math.sqrt(1 - cosalpha**2), cosalpha))
    beta = math.degrees(math.atan2(math.sqrt(1 - cosbeta**2), cosbeta))
    gamma = math.degrees(math.atan2(math.sqrt(1 - cosgamma**2), cosgamma))
    volume = abs(
        ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)
    )
    return a, b, c, alpha, beta, gamma, volume


def _ofmt(value: float) -> float:
    """Round through awk's default OFMT (%.6g), as the shell pipeline does."""
    return float(f"{value:.6g}")


def cell_row(cell, system: str) -> list[float]:
    """The CELL row for one configuration, in the columns of *system*."""
    key = normalise_system(system)
    # Rounded first: the shell scripts compute these columns from the %.6g
    # values in cellparams.dat, not from full-precision intermediates.
    _a, b, c, alpha, beta, gamma, v = (
        _ofmt(value) for value in cell_parameters(cell)
    )
    rad = math.radians
    if key == "cub":
        return [v ** (1 / 3), v]
    if key == "tet":
        return [math.sqrt(v / c), c, v]
    if key == "ort":
        return [v / (b * c), b, c, v]
    if key == "hex":
        # 0.866 (not sqrt(3)/2 to full precision) matches sod_vasp_cell.sh.
        return [math.sqrt(v / (0.866 * c)), c, v]
    if key == "rho":
        ca = math.cos(rad(alpha))
        return [(v / math.sqrt(1 - 3 * ca**2 + 2 * ca**3)) ** (1 / 3), alpha, v]
    if key == "mon":
        return [v / (b * c * math.sin(rad(beta))), b, c, beta, v]
    ca, cb, cg = math.cos(rad(alpha)), math.cos(rad(beta)), math.cos(rad(gamma))
    root = math.sqrt(1 - ca**2 - cb**2 - cg**2 + 2 * ca * cb * cg)
    return [v / (b * c * root), b, c, alpha, beta, gamma, v]


def format_row(values: list[float]) -> str:
    """Match awk's default %.6g output, which is what the shell scripts emit."""
    return " ".join(f"{value:.6g}" for value in values)


def write_cell_file(path, cells, system: str) -> None:
    """Write nXX/CELL: one row per configuration, in cYY order, no index column.

    Rows are positional -- like DATA, and unlike ENERGIES there is no
    configuration index -- so *cells* must cover every configuration of the
    level, in order.
    """
    key = normalise_system(system)
    lines = [format_row(cell_row(cell, key)) + "\n" for cell in cells]
    path.write_text("".join(lines))
