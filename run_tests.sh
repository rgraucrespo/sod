#!/usr/bin/env bash
# run_tests.sh — SOD regression test suite
#
# Usage: ./run_tests.sh
#
# Each test runs in a fresh temp directory so committed reference files are
# never overwritten.  Exit status is 0 if all tests pass, 1 otherwise.
#
# combsod tests (example02-14):
#   Copy INSOD + SGO to tmpdir, run combsod, diff every n*/OUTSOD against
#   the committed reference.
#
# genersod tests (example01/FILER*):
#   Copy input + template files to tmpdir, run combsod then genersod, check
#   for errors, and diff the first generated structure file against the
#   committed reference.
#
# statsod tests (canonical statistics):
#   Copy a single n*/OUTSOD + ENERGIES + DATA + SPECTRA + XSPEC to tmpdir,
#   run statsod, diff thermodynamics.dat, ave_data.dat, and ave_spectra.dat
#   against committed references.
#
# gcstatsod tests (grand-canonical statistics):
#   Copy the required range of n*/OUTSOD+ENERGIES+DATA+SPECTRA plus an
#   x????/ folder containing INGC to tmpdir, run sod_gcstat.sh, diff
#   thermodynamics.dat, ave_data.dat, and ave_spectra.dat.

ROOT="$(cd "$(dirname "$0")" && pwd)"
BIN="$ROOT/bin"
EX="$ROOT/examples"

pass=0; fail=0; skip=0

# ── helpers ──────────────────────────────────────────────────────────────────

label_fmt="%-42s"

pass_line()  { printf "PASS  $label_fmt\n" "$1"; }
fail_line()  { printf "FAIL  $label_fmt  %s\n" "$1" "$2"; }
skip_line()  { printf "SKIP  $label_fmt  %s\n" "$1" "$2"; }

indent() { sed 's/^/        /'; }

# ── test_combsod ─────────────────────────────────────────────────────────────
# $1 = display label   $2 = example directory
test_combsod() {
    local label="$1" dir="$2"

    # Require a non-empty INSOD
    if [ ! -s "$dir/INSOD" ]; then
        skip_line "$label" "(empty INSOD)"
        skip=$((skip+1)); return
    fi

    # Fresh working directory — only combsod inputs
    local tmp; tmp=$(mktemp -d)
    cp "$dir/INSOD" "$dir/SGO" "$tmp/"

    # Run combsod
    local out; out=$(cd "$tmp" && PATH="$BIN:$PATH" combsod 2>&1)
    local rc=$?
    if [ $rc -ne 0 ] || echo "$out" | grep -qi "error"; then
        fail_line "$label" "[combsod error]"
        echo "$out" | grep -i "error" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # For each OUTSOD combsod generated, compare against its reference.
    # (Direction: generated → reference, so stale committed files are ignored.)
    local gens mismatch=0
    gens=$(find "$tmp" -name OUTSOD 2>/dev/null | sort)
    if [ -z "$gens" ]; then
        fail_line "$label" "[no OUTSOD generated]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    while IFS= read -r gen; do
        local rel="${gen#$tmp/}"
        local ref="$dir/$rel"
        if [ ! -f "$ref" ]; then
            fail_line "$label" "[$rel has no reference]"
            mismatch=1
        elif ! diff -q "$ref" "$gen" >/dev/null 2>&1; then
            fail_line "$label" "[$rel differs]"
            diff "$ref" "$gen" | head -6 | indent
            mismatch=1
        fi
    done <<< "$gens"

    rm -rf "$tmp"
    if [ $mismatch -eq 0 ]; then
        pass_line "$label"; pass=$((pass+1))
    else
        fail=$((fail+1))
    fi
}

# ── test_genersod ─────────────────────────────────────────────────────────────
# $1 = display label   $2 = example directory   $3 = relative path of one
#     structure file to diff against its committed reference (may be empty)
test_genersod() {
    local label="$1" dir="$2" struct_rel="$3"

    if [ ! -s "$dir/INSOD" ]; then
        skip_line "$label" "(empty INSOD)"
        skip=$((skip+1)); return
    fi

    # Fresh working directory: copy everything except combsod/genersod outputs
    local tmp; tmp=$(mktemp -d)
    rsync -a \
        --exclude='n*/'        \
        --exclude='EQMATRIX'   \
        --exclude='filer'      \
        --exclude='OPERATORS'  \
        --exclude='supercell.cif' \
        --exclude='job_sender' \
        "$dir/" "$tmp/" 2>/dev/null

    # Run combsod (generates supercell.cif, filer, n*/OUTSOD, ...)
    local out; out=$(cd "$tmp" && PATH="$BIN:$PATH" combsod 2>&1)
    local rc=$?
    if [ $rc -ne 0 ] || echo "$out" | grep -qi "error"; then
        fail_line "$label" "[combsod error]"
        echo "$out" | grep -i "error" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Run genersod
    out=$(cd "$tmp" && PATH="$BIN:$PATH" genersod 2>&1)
    rc=$?
    if [ $rc -ne 0 ] || echo "$out" | grep -qi "error"; then
        fail_line "$label" "[genersod error]"
        echo "$out" | grep -i "error" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Diff the nominated structure file against the committed reference
    if [ -n "$struct_rel" ] && [ -f "$dir/$struct_rel" ]; then
        if ! diff -q "$dir/$struct_rel" "$tmp/$struct_rel" >/dev/null 2>&1; then
            fail_line "$label" "[$struct_rel differs]"
            diff "$dir/$struct_rel" "$tmp/$struct_rel" | head -8 | indent
            fail=$((fail+1)); rm -rf "$tmp"; return
        fi
    fi

    pass_line "$label"; pass=$((pass+1))
    rm -rf "$tmp"
}

# ── test_stat ────────────────────────────────────────────────────────────────
# $1 = display label   $2 = n-level directory (must contain OUTSOD, ENERGIES;
#      optionally DATA and SPECTRA)
# Reference files that must already be committed in $2:
#   thermodynamics.dat   (always)
#   ave_data.dat         (if DATA is present)
#   ave_spectra.dat      (if SPECTRA is present)
test_stat() {
    local label="$1" ndir="$2"

    if [ ! -f "$ndir/OUTSOD" ] || [ ! -f "$ndir/ENERGIES" ]; then
        skip_line "$label" "(missing OUTSOD or ENERGIES)"
        skip=$((skip+1)); return
    fi

    local tmp; tmp=$(mktemp -d)
    cp "$ndir/OUTSOD" "$ndir/ENERGIES" "$tmp/"
    [ -f "$ndir/DATA"    ] && cp "$ndir/DATA"    "$tmp/"
    [ -f "$ndir/SPECTRA" ] && cp "$ndir/SPECTRA" "$tmp/"
    [ -f "$ndir/XSPEC"   ] && cp "$ndir/XSPEC"   "$tmp/"

    local out; out=$(cd "$tmp" && PATH="$BIN:$PATH" statsod 2>&1)
    local rc=$?
    if [ $rc -ne 0 ] || echo "$out" | grep -qi "error"; then
        fail_line "$label" "[statsod error]"
        echo "$out" | grep -i "error" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local mismatch=0
    for f in thermodynamics.dat ave_data.dat ave_spectra.dat; do
        [ ! -f "$tmp/$f" ] && continue          # not produced (no DATA/SPECTRA)
        local ref="$ndir/$f"
        if [ ! -f "$ref" ]; then
            fail_line "$label" "[$f has no reference]"
            mismatch=1
        elif ! diff -q "$ref" "$tmp/$f" >/dev/null 2>&1; then
            fail_line "$label" "[$f differs]"
            diff "$ref" "$tmp/$f" | head -6 | indent
            mismatch=1
        fi
    done

    rm -rf "$tmp"
    if [ $mismatch -eq 0 ]; then
        pass_line "$label"; pass=$((pass+1))
    else
        fail=$((fail+1))
    fi
}

# ── test_gcstat ───────────────────────────────────────────────────────────────
# $1 = display label
# $2 = main example directory (contains the n*/ subdirectories)
# $3 = name of the x????/ subfolder within $2 (contains INGC and references)
# Reference files that must already be committed in $2/$3:
#   thermodynamics.dat, ave_data.dat, ave_spectra.dat
test_gcstat() {
    local label="$1" maindir="$2" xname="$3"
    local xdir="$maindir/$xname"

    if [ ! -f "$xdir/INGC" ]; then
        skip_line "$label" "(no INGC in $xname)"
        skip=$((skip+1)); return
    fi

    # Read the substitution range from INGC (line 2, skip comment)
    local nsubsmin nsubsmax
    read nsubsmin nsubsmax < <(grep -v "^#" "$xdir/INGC" | head -1)

    local tmp; tmp=$(mktemp -d)

    # Copy n*/OUTSOD+ENERGIES+DATA+SPECTRA for the required range
    local missing=0
    for ((i=nsubsmin; i<=nsubsmax; i++)); do
        local tag; tag=$(printf "%02d" $i)
        local nsrc="$maindir/n${tag}"
        if [ ! -f "$nsrc/OUTSOD" ] || [ ! -f "$nsrc/ENERGIES" ]; then
            fail_line "$label" "[n${tag}/OUTSOD or ENERGIES missing]"
            missing=1; break
        fi
        mkdir -p "$tmp/n${tag}"
        cp "$nsrc/OUTSOD" "$nsrc/ENERGIES" "$tmp/n${tag}/"
        [ -f "$nsrc/DATA"    ] && cp "$nsrc/DATA"    "$tmp/n${tag}/"
        [ -f "$nsrc/SPECTRA" ] && cp "$nsrc/SPECTRA" "$tmp/n${tag}/"
    done
    if [ $missing -ne 0 ]; then
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Copy the x????/ folder (INGC + optional TEMPERATURES/XSPEC)
    mkdir -p "$tmp/$xname"
    cp "$xdir/INGC" "$tmp/$xname/"
    [ -f "$xdir/TEMPERATURES" ] && cp "$xdir/TEMPERATURES" "$tmp/$xname/"
    [ -f "$xdir/XSPEC"        ] && cp "$xdir/XSPEC"        "$tmp/$xname/"

    local out; out=$(cd "$tmp/$xname" && PATH="$BIN:$PATH" sod_gcstat.sh 2>&1)
    local rc=$?
    if [ $rc -ne 0 ] || echo "$out" | grep -qi "error"; then
        fail_line "$label" "[sod_gcstat.sh error]"
        echo "$out" | grep -i "error" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local mismatch=0
    for f in thermodynamics.dat ave_data.dat ave_spectra.dat; do
        [ ! -f "$tmp/$xname/$f" ] && continue
        local ref="$xdir/$f"
        if [ ! -f "$ref" ]; then
            fail_line "$label" "[$f has no reference]"
            mismatch=1
        elif ! diff -q "$ref" "$tmp/$xname/$f" >/dev/null 2>&1; then
            fail_line "$label" "[$f differs]"
            diff "$ref" "$tmp/$xname/$f" | head -6 | indent
            mismatch=1
        fi
    done

    rm -rf "$tmp"
    if [ $mismatch -eq 0 ]; then
        pass_line "$label"; pass=$((pass+1))
    else
        fail=$((fail+1))
    fi
}

# ── combsod tests ─────────────────────────────────────────────────────────────

echo "SOD regression tests"
echo "===================="
echo ""
printf "combsod  (examples 02-14)\n"
printf "%s\n" "-------------------------------------------"

for ex in "$EX"/example??; do
    [ "$(basename "$ex")" = "example01" ] && continue  # handled in genersod section
    test_combsod "$(basename "$ex")" "$ex"
done

# ── genersod tests ────────────────────────────────────────────────────────────

echo ""
printf "genersod  (example01 FILER subdirectories)\n"
printf "%s\n" "-------------------------------------------"

test_genersod "example01/FILER1_gulp"    "$EX/example01/FILER1_gulp"    "n04/c01/input.gin"
test_genersod "example01/FILER2_lammps"  "$EX/example01/FILER2_lammps"  "n04/c01/conf.data"
test_genersod "example01/FILER11_vasp"   "$EX/example01/FILER11_vasp"   "n04/c01/POSCAR"
test_genersod "example01/FILER12_castep" "$EX/example01/FILER12_castep" "n04/c01/castep.cell"
test_genersod "example01/FILER13_QE"     "$EX/example01/FILER13_QE"     "n04/c01/pw.in"

# ── statsod tests (canonical) ────────────────────────────────────────────────

echo ""
printf "statsod  (canonical statistics)\n"
printf "%s\n" "-------------------------------------------"

test_stat "example05/n02 (energy+DATA+SPECTRA)" "$EX/example05/n02"

# ── gcstatsod tests (grand-canonical) ────────────────────────────────────────

echo ""
printf "gcstatsod  (grand-canonical statistics)\n"
printf "%s\n" "-------------------------------------------"

test_gcstat "example05/test_gcstat (n10-n16, x=0.875)" "$EX/example05" "test_gcstat"

# ── summary ───────────────────────────────────────────────────────────────────

echo ""
echo "Results: $pass passed, $fail failed, $skip skipped"
[ $fail -eq 0 ]
