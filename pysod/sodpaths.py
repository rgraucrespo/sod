"""SOD directory conventions.

Python port of the discovery helpers in ``bin/sod_common.sh``.  The rules must
stay in step with that file: a SODPROJECT is a directory holding both INSOD and
SGO, level directories are ``nXX`` (``n04``, ``n02_02``, ``n001_001``), and
configuration directories are ``cYY`` with a *dynamic* zero padding chosen by
genersod from the configuration count.
"""

from __future__ import annotations

import re
from pathlib import Path

# Mirrors sod_find_enclosing_level_name in bin/sod_common.sh.
LEVEL_RE = re.compile(r"^n[0-9]+(_[0-9]+)*$")

# Per-configuration structure file written by genersod for each FILER value.
# Only the ASE-readable ones are usable as MLIP input.
FILER_STRUCTURE_FILES = {
    0: "configuration.cif",
    11: "POSCAR",
}

# ASE format name for each supported structure file, so writers do not have to
# guess from a suffix (POSCAR has none).
FILER_ASE_FORMAT = {
    0: "cif",
    11: "vasp",
}

FILER_DESCRIPTIONS = {
    -1: "no calculation files",
    0: "CIF",
    1: "GULP",
    2: "LAMMPS",
    11: "VASP",
    12: "CASTEP",
    13: "QE",
}


class SodLayoutError(Exception):
    """Raised when the on-disk layout does not match SOD's conventions."""


def find_project_root(start: Path | str = ".") -> Path:
    """Walk upward until a directory contains both INSOD and SGO."""
    directory = Path(start).resolve()
    for candidate in (directory, *directory.parents):
        if (candidate / "INSOD").is_file() and (candidate / "SGO").is_file():
            return candidate
    raise SodLayoutError(
        f"could not locate SODPROJECT/ from {directory}.\n"
        "       Searched upward for a directory containing INSOD and SGO."
    )


def find_level_name(root: Path, start: Path | str = ".") -> str | None:
    """Return the enclosing ``nXX`` directory name, or None if at the root."""
    root = Path(root).resolve()
    directory = Path(start).resolve()
    while directory != root and directory != directory.parent:
        if LEVEL_RE.match(directory.name):
            return directory.name
        directory = directory.parent
    return None


def level_dirs(root: Path) -> list[Path]:
    """All ``nXX`` level directories directly under a SODPROJECT, sorted."""
    return sorted(
        (path for path in root.glob("n*") if path.is_dir() and LEVEL_RE.match(path.name)),
        key=lambda path: path.name,
    )


def config_dirs(level_dir: Path) -> list[tuple[int, Path]]:
    """Return ``(m, directory)`` pairs for every ``cYY`` configuration directory.

    The zero padding of ``cYY`` depends on the configuration count (71 configs
    give ``c01``, 8043 give ``c0001``), so the index is recovered by a base-10
    integer parse -- the Python equivalent of ``$((10#${name#c}))``.  Entries
    such as the generated ``cSGO`` file are skipped by the digit test.
    """
    pairs = [
        (int(path.name[1:], 10), path)
        for path in level_dir.glob("c*")
        if path.is_dir() and path.name[1:].isdigit()
    ]
    return sorted(pairs, key=lambda pair: pair[0])


def ensemble_config_count(path: Path) -> int:
    """Return the configuration count from the first line of a v3 ENSEMBLE.

    Line 1 is ``<type> ensemble[ (T K)]: N configurations; ...``, so the count
    is the integer token immediately before ``configurations``.  Version 2
    files are not accepted -- their count is a bare integer on a later line,
    which cannot be told apart from a data row.
    """
    path = Path(path)
    if not path.is_file():
        raise SodLayoutError(f"ENSEMBLE file not found: {path}")
    header = next(
        (line for line in path.read_text().splitlines() if line.strip()), ""
    )
    words = header.split()
    if words and not words[0].startswith("#"):
        for position, word in enumerate(words[1:], start=1):
            if word.rstrip(";:") == "configurations" and words[position - 1].isdigit():
                count = int(words[position - 1])
                if count >= 1:
                    return count
    raise SodLayoutError(
        f"{path} is not a version 3 ENSEMBLE. Expected a first line such as "
        "'Enumerated ensemble: 71 configurations; sum_degeneracies = 4096'. "
        "Version 2 ENSEMBLE files are no longer supported; regenerate the level "
        "with sod_comb.sh or sod_random.sh"
    )


def read_filer(insod_path: Path) -> int:
    """FILER is the last non-blank, non-comment line of INSOD."""
    insod_path = Path(insod_path)
    if not insod_path.is_file():
        raise SodLayoutError(f"INSOD file not found: {insod_path}")
    last = ""
    for raw in insod_path.read_text().splitlines():
        line = raw.strip()
        if line and not line.startswith("#"):
            last = line
    if not last:
        raise SodLayoutError(f"could not locate FILER line in {insod_path}")
    try:
        return int(last)
    except ValueError as exc:
        raise SodLayoutError(
            f"could not parse FILER from INSOD tail line: {last}"
        ) from exc


def structure_filename(filer: int) -> str:
    """Map a FILER value to the per-configuration structure file to read."""
    try:
        return FILER_STRUCTURE_FILES[filer]
    except KeyError:
        description = FILER_DESCRIPTIONS.get(filer, "unknown")
        supported = ", ".join(
            f"{value} ({FILER_DESCRIPTIONS[value]})" for value in sorted(FILER_STRUCTURE_FILES)
        )
        raise SodLayoutError(
            f"FILER {filer} ({description}) does not produce an ASE-readable structure file.\n"
            f"       Supported: {supported}.\n"
            "       Use -structure NAME to point at a structure file explicitly."
        ) from None
