#!/usr/bin/env bash
# sod_run_tests.sh — SOD regression test suite
#
# Usage: ./bin/sod_run_tests.sh
#
# Each test runs in a fresh temp directory so committed reference files are
# never overwritten.  Exit status is 0 if all tests pass, 1 otherwise.
#
# combsod tests (example02-14):
#   Copy INSOD + SGO to tmpdir, run combsod, diff every n*/ENSEMBLE against
#   the committed reference.
#
# genersod tests (example01/FILER*):
#   Copy input + template files to tmpdir, run combsod then genersod, check
#   for errors, and diff the first generated structure file against the
#   committed reference.
#
# statsod tests (canonical statistics):
#   Copy a single n*/ENSEMBLE + ENERGIES + DATA + SPECTRA + XSPEC to tmpdir,
#   run statsod, diff thermodynamics.dat, ave_data.dat, and ave_spectra.dat
#   against committed references.
#
# gcstatsod tests (grand-canonical statistics):
#   Copy the required range of n*/ENSEMBLE+ENERGIES+DATA+SPECTRA plus an
#   x????/ folder containing INGC to tmpdir, run sod_gcstat.sh, diff
#   thermodynamics.dat, ave_data.dat, and ave_spectra.dat.

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
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
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[combsod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # For each ENSEMBLE combsod generated, compare against its reference.
    # (Direction: generated → reference, so stale committed files are ignored.)
    local gens mismatch=0
    gens=$(find "$tmp" -name ENSEMBLE 2>/dev/null | sort)
    if [ -z "$gens" ]; then
        fail_line "$label" "[no ENSEMBLE generated]"
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

    # Run combsod (generates supercell.cif, filer, n*/ENSEMBLE, ...)
    local out; out=$(cd "$tmp" && PATH="$BIN:$PATH" combsod 2>&1)
    local rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[combsod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Run genersod
    out=$(cd "$tmp" && PATH="$BIN:$PATH" genersod 2>&1)
    rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[genersod error]"
        echo "$out" | head -3 | indent
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
# $1 = display label   $2 = n-level directory (must contain ENSEMBLE, ENERGIES;
#      optionally TEMPERATURES, DATA and SPECTRA)
# If TEMPERATURES is present it is used; otherwise statsod's default 0/300/1000 K.
# Reference files that must already be committed in $2:
#   thermodynamics.dat   (always)
#   ave_data.dat         (if DATA is present)
#   ave_spectra.dat      (if SPECTRA is present)
test_stat() {
    local label="$1" ndir="$2"

    if [ ! -f "$ndir/ENSEMBLE" ] || [ ! -f "$ndir/ENERGIES" ]; then
        skip_line "$label" "(missing ENSEMBLE or ENERGIES)"
        skip=$((skip+1)); return
    fi

    local tmp; tmp=$(mktemp -d)
    cp "$ndir/ENSEMBLE" "$ndir/ENERGIES" "$tmp/"
    [ -f "$ndir/TEMPERATURES" ] && cp "$ndir/TEMPERATURES" "$tmp/"
    [ -f "$ndir/DATA"    ] && cp "$ndir/DATA"    "$tmp/"
    [ -f "$ndir/SPECTRA" ] && cp "$ndir/SPECTRA" "$tmp/"
    [ -f "$ndir/XSPEC"   ] && cp "$ndir/XSPEC"   "$tmp/"

    local out; out=$(cd "$tmp" && PATH="$BIN:$PATH" statsod 2>&1)
    local rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[statsod error]"
        echo "$out" | head -3 | indent
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

# ── test_pmesod ───────────────────────────────────────────────────────────────
# $1 = display label
# $2 = main example directory (contains INSOD, SGO, n00-n03, n21-n24)
# $3 = target level
# Reference file must already be committed in $2/pme_test_ref/ENERGIES
test_pmesod() {
    local label="$1" maindir="$2" tlvl="$3"

    if [ ! -f "$maindir/INSOD" ] || [ ! -f "$maindir/SGO" ]; then
        skip_line "$label" "(missing INSOD or SGO)"
        skip=$((skip+1)); return
    fi

    # Low-side levels: 0-3; high-side levels: 21-24

    # Check that all required level directories and data files exist
    local missing=0
    for i in 0 1 2 3 21 22 23 24; do
        local tag; tag=$(printf "%02d" $i)
        local ndir="$maindir/n${tag}"
        if [ ! -f "$ndir/ENSEMBLE" ] || [ ! -f "$ndir/ENERGIES" ]; then
            fail_line "$label" "[n${tag}/ENSEMBLE or ENERGIES missing]"
            missing=1; break
        fi
    done
    if [ $missing -ne 0 ]; then
        fail=$((fail+1)); return
    fi

    local tmp; tmp=$(mktemp -d)
    cp "$maindir/INSOD" "$maindir/SGO" "$tmp/"

    # Copy low and high level reference data
    for i in 0 1 2 3 21 22 23 24; do
        local tag; tag=$(printf "%02d" $i)
        mkdir -p "$tmp/n${tag}"
        cp "$maindir/n${tag}/ENSEMBLE" "$maindir/n${tag}/ENERGIES" "$tmp/n${tag}/"
    done

    # Run combsod to generate EQMATRIX
    local out; out=$(cd "$tmp" && PATH="$BIN:$PATH" combsod 2>&1)
    local rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[combsod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Run pmesod
    out=$(cd "$tmp" && PATH="$BIN:$PATH" pmesod 2>&1)
    rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[pmesod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Check PMEh ENERGIES from target level
    local ttag; ttag=$(printf "%02d" $tlvl)
    local ref="$maindir/pme_test_ref/ENERGIES"
    local gen="$tmp/n${ttag}/PMEh/ENERGIES"
    if [ ! -f "$ref" ]; then
        fail_line "$label" "[ENERGIES has no reference]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    elif ! diff -q "$ref" "$gen" >/dev/null 2>&1; then
        fail_line "$label" "[ENERGIES differs]"
        diff "$ref" "$gen" | head -6 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    rm -rf "$tmp"
    pass_line "$label"; pass=$((pass+1))
}

# ── test_pmesod_example17 ─────────────────────────────────────────────────────
# example17: n00–n03 (low, order_cap=3) + n24–n27 (high), target n04/PMEh.
# pme.model is copied so pmesod sees pme_ch=2/order_cap=3/n_calib=0/alpha=2/eta=0.
# n04 is the prediction target, not a training level.
# n04/ENSEMBLE and n04/ENERGIES are deliberately included in the tmp dir to
# verify that pmesod does NOT use them as training data (the target-level guard).
# Without the guard, pmesod would crash with a "training data conflict" error.
# n_calib=0 → manual eps mode (eps=1, alpha=2, eta=0): "epsilon read/fitted from pme.model."
# $1 = display label
# $2 = main example directory
test_pmesod_example17() {
    local label="$1" maindir="$2"

    if [ ! -f "$maindir/INSOD" ] || [ ! -f "$maindir/SGO" ] || [ ! -f "$maindir/pme.model" ]; then
        skip_line "$label" "(missing INSOD, SGO or pme.model)"
        skip=$((skip+1)); return
    fi

    local missing=0
    for i in 0 1 2 3 4 24 25 26 27; do
        local tag; tag=$(printf "%02d" $i)
        if [ ! -f "$maindir/n${tag}/ENSEMBLE" ] || [ ! -f "$maindir/n${tag}/ENERGIES" ]; then
            fail_line "$label" "[n${tag}/ENSEMBLE or ENERGIES missing]"
            missing=1; break
        fi
    done
    [ $missing -ne 0 ] && { fail=$((fail+1)); return; }

    local tmp; tmp=$(mktemp -d)
    cp "$maindir/INSOD" "$maindir/SGO" "$maindir/pme.model" "$tmp/"

    # Include n04 so the target-level guard is exercised: pmesod must not use
    # n04/ENERGIES as a training reference even though the file is present.
    for i in 0 1 2 3 4 24 25 26 27; do
        local tag; tag=$(printf "%02d" $i)
        mkdir -p "$tmp/n${tag}"
        cp "$maindir/n${tag}/ENSEMBLE" "$maindir/n${tag}/ENERGIES" "$tmp/n${tag}/"
    done

    local out; out=$(cd "$tmp" && PATH="$BIN:$PATH" combsod 2>&1)
    local rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[combsod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    out=$(cd "$tmp" && PATH="$BIN:$PATH" pmesod 2>&1)
    rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[pmesod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local ref="$maindir/pme_test_ref/ENERGIES"
    local gen="$tmp/n04/PMEh/ENERGIES"
    if [ ! -f "$ref" ]; then
        fail_line "$label" "[ENERGIES has no reference]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    elif ! diff -q "$ref" "$gen" >/dev/null 2>&1; then
        fail_line "$label" "[ENERGIES differs]"
        diff "$ref" "$gen" | head -6 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    rm -rf "$tmp"
    pass_line "$label"; pass=$((pass+1))
}

# ── test_mcsod ────────────────────────────────────────────────────────────────
# $1 = display label
# $2 = main example directory (contains INSOD, SGO, INMC, n00-n03, n21-n24)
# $3 = target level (where OUTMC is written)
# Reference file must already be committed in $2/pme_ref/OUTMC
test_mcsod() {
    local label="$1" maindir="$2" tlvl="$3"

    if [ ! -f "$maindir/INSOD" ] || [ ! -f "$maindir/SGO" ] || [ ! -f "$maindir/INMC" ]; then
        skip_line "$label" "(missing INSOD, SGO, or INMC)"
        skip=$((skip+1)); return
    fi

    # Low-side levels: 0-3; high-side levels: 21-24

    # Check that all required level directories and data files exist
    local missing=0
    for i in 0 1 2 3 21 22 23 24; do
        local tag; tag=$(printf "%02d" $i)
        local ndir="$maindir/n${tag}"
        if [ ! -f "$ndir/ENSEMBLE" ] || [ ! -f "$ndir/ENERGIES" ]; then
            fail_line "$label" "[n${tag}/ENSEMBLE or ENERGIES missing]"
            missing=1; break
        fi
    done
    if [ $missing -ne 0 ]; then
        fail=$((fail+1)); return
    fi

    local tmp; tmp=$(mktemp -d)
    cp "$maindir/INSOD" "$maindir/SGO" "$maindir/INMC" "$tmp/"
    printf '300.0\n' > "$tmp/TEMPERATURES"

    # Copy low and high level reference data
    for i in 0 1 2 3 21 22 23 24; do
        local tag; tag=$(printf "%02d" $i)
        mkdir -p "$tmp/n${tag}"
        cp "$maindir/n${tag}/ENSEMBLE" "$maindir/n${tag}/ENERGIES" "$tmp/n${tag}/"
    done

    # Run combsod then pmesod to generate Hamiltonian (needed for mcsod)
    local out; out=$(cd "$tmp" && PATH="$BIN:$PATH" combsod 2>&1)
    local rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[combsod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    out=$(cd "$tmp" && PATH="$BIN:$PATH" pmesod 2>&1)
    rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[pmesod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Run MC through the script so TEMPERATURES is looped outside mcsod
    out=$(cd "$tmp" && PATH="$BIN:$PATH" sod_mc.sh 2>&1)
    rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[mcsod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Check OUTMC from target level
    local ttag; ttag=$(printf "%02d" $tlvl)
    local ref="$maindir/pme_test_ref/OUTMC"
    local gen="$tmp/n${ttag}/PMEh/MCT_300K/OUTMC"
    if [ ! -f "$ref" ]; then
        fail_line "$label" "[OUTMC has no reference]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    elif ! diff -q "$ref" "$gen" >/dev/null 2>&1; then
        fail_line "$label" "[OUTMC differs]"
        diff "$ref" "$gen" | head -6 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    rm -rf "$tmp"
    pass_line "$label"; pass=$((pass+1))
}

# ── test_mcstat ───────────────────────────────────────────────────────────────
# $1 = display label
# $2 = main example directory (contains INSOD, SGO, INMC, n00-n03, n21-n24)
# $3 = target level (where the MCT_*K/ directories are written)
# Runs the MC chain over a multi-temperature ladder, then thermodynamic
# integration (mcstatsod) at the target level.
# Reference file must already be committed in $2/pme_test_ref/thermodynamics.dat
test_mcstat() {
    local label="$1" maindir="$2" tlvl="$3"

    if [ ! -f "$maindir/INSOD" ] || [ ! -f "$maindir/SGO" ] || [ ! -f "$maindir/INMC" ]; then
        skip_line "$label" "(missing INSOD, SGO, or INMC)"
        skip=$((skip+1)); return
    fi

    # Low-side levels: 0-3; high-side levels: 21-24
    local missing=0
    for i in 0 1 2 3 21 22 23 24; do
        local tag; tag=$(printf "%02d" $i)
        local ndir="$maindir/n${tag}"
        if [ ! -f "$ndir/ENSEMBLE" ] || [ ! -f "$ndir/ENERGIES" ]; then
            fail_line "$label" "[n${tag}/ENSEMBLE or ENERGIES missing]"
            missing=1; break
        fi
    done
    if [ $missing -ne 0 ]; then
        fail=$((fail+1)); return
    fi

    local tmp; tmp=$(mktemp -d)
    cp "$maindir/INSOD" "$maindir/SGO" "$maindir/INMC" "$tmp/"
    # Multi-temperature ladder — thermodynamic integration needs >= 2 temperatures
    printf '1000.0\n600.0\n300.0\n' > "$tmp/TEMPERATURES"

    # Copy low and high level reference data
    for i in 0 1 2 3 21 22 23 24; do
        local tag; tag=$(printf "%02d" $i)
        mkdir -p "$tmp/n${tag}"
        cp "$maindir/n${tag}/ENSEMBLE" "$maindir/n${tag}/ENERGIES" "$tmp/n${tag}/"
    done

    # Run combsod then pmesod to generate the Hamiltonian (needed for mcsod)
    local out; out=$(cd "$tmp" && PATH="$BIN:$PATH" combsod 2>&1)
    local rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[combsod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    out=$(cd "$tmp" && PATH="$BIN:$PATH" pmesod 2>&1)
    rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[pmesod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Run MC through the script so TEMPERATURES is looped outside mcsod
    out=$(cd "$tmp" && PATH="$BIN:$PATH" sod_mc.sh 2>&1)
    rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[mcsod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Thermodynamic integration over the sampled temperatures
    local ttag; ttag=$(printf "%02d" $tlvl)
    local pmedir="$tmp/n${ttag}/PMEh"
    if [ ! -d "$pmedir" ]; then
        fail_line "$label" "[n${ttag}/PMEh not generated]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    out=$(cd "$pmedir" && PATH="$BIN:$PATH" sod_mcstat.sh 2>&1)
    rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[mcstatsod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Check thermodynamics.dat against the committed reference
    local ref="$maindir/pme_test_ref/thermodynamics.dat"
    local gen="$pmedir/thermodynamics.dat"
    if [ ! -f "$ref" ]; then
        fail_line "$label" "[thermodynamics.dat has no reference]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    elif ! diff -q "$ref" "$gen" >/dev/null 2>&1; then
        fail_line "$label" "[thermodynamics.dat differs]"
        diff "$ref" "$gen" | head -6 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    rm -rf "$tmp"
    pass_line "$label"; pass=$((pass+1))
}

# ── test_mcstat_vs_enum ───────────────────────────────────────────────────────
# Cross-check: thermodynamic integration (mcstatsod) vs the exact full
# enumeration (statsod) on the *same* PMEh Hamiltonian. Both paths use the PMEh
# energies that pmesod predicts for every configuration of the target level, so
# any disagreement is pure method error (MC sampling + TI discretization).
# Asserts max|F_TI(T) - F_exact(T)| <= tol over the temperature ladder.
# $1 = display label  $2 = main example dir  $3 = target level  $4 = tol (eV)
test_mcstat_vs_enum() {
    local label="$1" maindir="$2" tlvl="$3" tol="$4"

    if [ ! -f "$maindir/INSOD" ] || [ ! -f "$maindir/SGO" ] || [ ! -f "$maindir/INMC" ]; then
        skip_line "$label" "(missing INSOD, SGO, or INMC)"
        skip=$((skip+1)); return
    fi

    local missing=0
    for i in 0 1 2 3 21 22 23 24; do
        local tag; tag=$(printf "%02d" $i)
        if [ ! -f "$maindir/n${tag}/ENSEMBLE" ] || [ ! -f "$maindir/n${tag}/ENERGIES" ]; then
            fail_line "$label" "[n${tag}/ENSEMBLE or ENERGIES missing]"
            missing=1; break
        fi
    done
    if [ $missing -ne 0 ]; then fail=$((fail+1)); return; fi

    local tmp; tmp=$(mktemp -d)
    cp "$maindir/INSOD" "$maindir/SGO" "$maindir/INMC" "$tmp/"
    # Temperature ladder spanning the thermally interesting range. The lowest T
    # is the worst case for TI (integration edge); more points tighten agreement.
    printf '%s\n' 4000 2000 1200 800 500 300 > "$tmp/TEMPERATURES"
    for i in 0 1 2 3 21 22 23 24; do
        local tag; tag=$(printf "%02d" $i)
        mkdir -p "$tmp/n${tag}"
        cp "$maindir/n${tag}/ENSEMBLE" "$maindir/n${tag}/ENERGIES" "$tmp/n${tag}/"
    done

    # combsod (full enumeration of target level) + pmesod (PMEh energies)
    local out rc
    out=$(cd "$tmp" && PATH="$BIN:$PATH" combsod 2>&1); rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[combsod error]"; echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    out=$(cd "$tmp" && PATH="$BIN:$PATH" pmesod 2>&1); rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[pmesod error]"; echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local ttag; ttag=$(printf "%02d" $tlvl)
    local pmedir="$tmp/n${ttag}/PMEh"
    if [ ! -f "$pmedir/ENERGIES" ] || [ ! -f "$tmp/n${ttag}/ENSEMBLE" ]; then
        fail_line "$label" "[n${ttag} enumeration/ENERGIES not generated]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Exact path: statsod over the full enumeration with the PMEh energies.
    # pmesod writes single-column energies in configuration order; statsod wants
    # the two-column "m E" form, so prepend the 1-based index.
    mkdir -p "$tmp/enum"
    cp "$tmp/n${ttag}/ENSEMBLE" "$tmp/enum/ENSEMBLE"
    cp "$tmp/TEMPERATURES"      "$tmp/enum/TEMPERATURES"
    awk '{printf "%d  %s\n", NR, $1}' "$pmedir/ENERGIES" > "$tmp/enum/ENERGIES"
    out=$(cd "$tmp/enum" && PATH="$BIN:$PATH" statsod 2>&1); rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[statsod error]"; echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Approximate path: mcsod MC sampling (same PMEh) + mcstatsod TI
    out=$(cd "$tmp" && PATH="$BIN:$PATH" sod_mc.sh 2>&1); rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[mcsod error]"; echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    out=$(cd "$pmedir" && PATH="$BIN:$PATH" sod_mcstat.sh 2>&1); rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[mcstatsod error]"; echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Compare the F column at matching temperatures (keyed on integer T).
    local res maxd maxT npts verdict
    res=$(awk -v tol="$tol" '
        function ti(x) { return sprintf("%d", x + 0) }
        FNR == NR { if ($1 ~ /^[0-9]/) f_ex[ti($1)] = $3; next }
        ($1 ~ /^[0-9]/) {
            k = ti($1)
            if (!(k in f_ex)) next
            d = $3 - f_ex[k]; if (d < 0) d = -d
            if (d > maxd) { maxd = d; maxT = k }
            n++
        }
        END {
            if (n == 0) { print "0 0 0 ERR"; exit }
            printf "%.6f %s %d %s\n", maxd, maxT, n, (maxd <= tol ? "PASS" : "FAIL")
        }
    ' "$tmp/enum/thermodynamics.dat" "$pmedir/thermodynamics.dat")
    read -r maxd maxT npts verdict <<< "$res"

    rm -rf "$tmp"

    if [ "$verdict" = "PASS" ]; then
        pass_line "$label"
        printf "        max|dF| = %s eV at T=%s K over %s temps (tol %s eV)\n" \
            "$maxd" "$maxT" "$npts" "$tol"
        pass=$((pass+1))
    elif [ "$verdict" = "ERR" ]; then
        fail_line "$label" "[no matching temperatures between exact and TI output]"
        fail=$((fail+1))
    else
        fail_line "$label" "[max|dF|=$maxd eV at T=$maxT K over $npts temps exceeds tol $tol eV]"
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

    # Copy n*/ENSEMBLE+ENERGIES+DATA+SPECTRA for the required range
    local missing=0
    for ((i=nsubsmin; i<=nsubsmax; i++)); do
        local tag; tag=$(printf "%02d" $i)
        local nsrc="$maindir/n${tag}"
        if [ ! -f "$nsrc/ENSEMBLE" ] || [ ! -f "$nsrc/ENERGIES" ]; then
            fail_line "$label" "[n${tag}/ENSEMBLE or ENERGIES missing]"
            missing=1; break
        fi
        mkdir -p "$tmp/n${tag}"
        cp "$nsrc/ENSEMBLE" "$nsrc/ENERGIES" "$tmp/n${tag}/"
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
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[sod_gcstat.sh error]"
        echo "$out" | head -3 | indent
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

# ── test_sqssod ───────────────────────────────────────────────────────────────
# $1 = display label
# $2 = main example directory (contains INSOD, EQMATRIX, supercell.cif)
# $3 = nXX directory name (e.g. n08)
# Reference: $2/$3/OUTSQS (must be committed)
test_sqssod() {
    local label="$1" maindir="$2" nxx="$3"
    local ndir="$maindir/$nxx"

    for f in INSOD EQMATRIX supercell.cif; do
        if [ ! -f "$maindir/$f" ]; then
            skip_line "$label" "(missing $f)"
            skip=$((skip+1)); return
        fi
    done
    if [ ! -f "$ndir/ENSEMBLE" ]; then
        skip_line "$label" "(missing $nxx/ENSEMBLE)"
        skip=$((skip+1)); return
    fi
    if [ ! -f "$ndir/OUTSQS" ]; then
        fail_line "$label" "[$nxx/OUTSQS has no reference]"
        fail=$((fail+1)); return
    fi

    local tmp; tmp=$(mktemp -d)
    cp "$maindir/INSOD" "$maindir/EQMATRIX" "$maindir/supercell.cif" "$tmp/"
    [ -f "$maindir/INSQS" ] && cp "$maindir/INSQS" "$tmp/"
    mkdir -p "$tmp/$nxx"
    cp "$ndir/ENSEMBLE" "$tmp/$nxx/"
    [ -f "$ndir/INSQS" ] && cp "$ndir/INSQS" "$tmp/$nxx/"

    local out; out=$(cd "$tmp" && PATH="$BIN:$PATH" sqssod "$nxx" 2>&1)
    local rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[sqssod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local ref="$ndir/OUTSQS"
    local gen="$tmp/$nxx/OUTSQS"
    if ! diff -q "$ref" "$gen" >/dev/null 2>&1; then
        fail_line "$label" "[OUTSQS differs]"
        diff "$ref" "$gen" | head -6 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    rm -rf "$tmp"
    pass_line "$label"; pass=$((pass+1))
}

# ── test_gqssod ───────────────────────────────────────────────────────────────
# $1 = display label
# $2 = main example directory (contains INSOD, EQMATRIX, supercell.cif)
# $3 = nXX directory name (e.g. n08)
# Reference: $2/$3/OUTGQS (must be committed)
test_gqssod() {
    local label="$1" maindir="$2" nxx="$3"
    local ndir="$maindir/$nxx"

    for f in INSOD EQMATRIX supercell.cif; do
        if [ ! -f "$maindir/$f" ]; then
            skip_line "$label" "(missing $f)"
            skip=$((skip+1)); return
        fi
    done
    for f in ENSEMBLE ENERGIES OUTSQS; do
        if [ ! -f "$ndir/$f" ]; then
            skip_line "$label" "(missing $nxx/$f)"
            skip=$((skip+1)); return
        fi
    done
    if [ ! -f "$ndir/TEMPERATURES" ] && [ ! -f "$maindir/TEMPERATURES" ]; then
        skip_line "$label" "(missing TEMPERATURES)"
        skip=$((skip+1)); return
    fi
    if [ ! -f "$ndir/OUTGQS" ]; then
        fail_line "$label" "[$nxx/OUTGQS has no reference]"
        fail=$((fail+1)); return
    fi

    local tmp; tmp=$(mktemp -d)
    cp "$maindir/INSOD" "$maindir/EQMATRIX" "$maindir/supercell.cif" "$tmp/"
    [ -f "$maindir/INSQS" ] && cp "$maindir/INSQS" "$tmp/"
    [ -f "$maindir/TEMPERATURES" ] && cp "$maindir/TEMPERATURES" "$tmp/"
    mkdir -p "$tmp/$nxx"
    cp "$ndir/ENSEMBLE" "$ndir/ENERGIES" "$ndir/OUTSQS" "$tmp/$nxx/"
    [ -f "$ndir/INSQS" ] && cp "$ndir/INSQS" "$tmp/$nxx/"
    [ -f "$ndir/TEMPERATURES" ] && cp "$ndir/TEMPERATURES" "$tmp/$nxx/"

    local out; out=$(cd "$tmp" && PATH="$BIN:$PATH" gqssod "$nxx" 2>&1)
    local rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[gqssod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local ref="$ndir/OUTGQS"
    local gen="$tmp/$nxx/OUTGQS"
    if ! diff -q "$ref" "$gen" >/dev/null 2>&1; then
        fail_line "$label" "[OUTGQS differs]"
        diff "$ref" "$gen" | head -6 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    rm -rf "$tmp"
    pass_line "$label"; pass=$((pass+1))
}

# ── combsod tests ─────────────────────────────────────────────────────────────

echo "SOD regression tests"
echo "===================="
echo ""
printf "combsod  (examples 02-14)\n"
printf "%s\n" "-------------------------------------------"

for ex in "$EX"/example0[2-9] "$EX"/example1[0-4]; do
    [ -d "$ex" ] && test_combsod "$(basename "$ex")" "$ex"
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
test_stat "example16/n08 (8043 configs, T=0/1000/1e6 K)" "$EX/example16/n08"

# ── pmesod tests (Periodic Motif Expansion) ────────────────────────────────────

echo ""
printf "pmesod  (Periodic Motif Expansion)\n"
printf "%s\n" "-------------------------------------------"

test_pmesod "example15/pmesod (n00-n03 low, n21-n24 high, target=n12)" "$EX/example15" 12
test_pmesod_example17 "example17/pmesod (n00-n03 low, n24-n27 high, target=n04, PMEh)" "$EX/example17"

# ── mcsod tests (Monte Carlo sampling) ──────────────────────────────────────────

echo ""
printf "mcsod  (Monte Carlo sampling)\n"
printf "%s\n" "-------------------------------------------"

test_mcsod "example15/mcsod (n12, 300K, 5000 samples)" "$EX/example15" 12

# ── mcstatsod tests (thermodynamic integration) ──────────────────────────────

echo ""
printf "mcstatsod  (thermodynamic integration)\n"
printf "%s\n" "-------------------------------------------"

test_mcstat "example15/mcstat (n12, T=300/600/1000K)" "$EX/example15" 12
test_mcstat_vs_enum "example15/mcstat vs full enum (F within tol)" "$EX/example15" 12 0.005

# ── gcstatsod tests (grand-canonical) ────────────────────────────────────────

echo ""
printf "gcstatsod  (grand-canonical statistics)\n"
printf "%s\n" "-------------------------------------------"

test_gcstat "example05/test_gcstat (n10-n16, x=0.875)" "$EX/example05" "test_gcstat"

# ── sqssod tests (Special Quasirandom Structures) ────────────────────────────

echo ""
printf "sqssod  (Special Quasirandom Structures)\n"
printf "%s\n" "-------------------------------------------"

test_sqssod "example16/sqssod (n08, 8043 configs)" "$EX/example16" "n08"

# ── gqssod tests (Generalized Quasirandom Structures) ────────────────────────

echo ""
printf "gqssod  (Generalized Quasirandom Structures)\n"
printf "%s\n" "-------------------------------------------"

test_gqssod "example16/gqssod (n08, T=0/1000/1e6 K)" "$EX/example16" "n08"

# ── summary ───────────────────────────────────────────────────────────────────

echo ""
echo "Results: $pass passed, $fail failed, $skip skipped"
[ $fail -eq 0 ]
