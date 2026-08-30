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

. "${BIN}/sod_common.sh"

pass=0; fail=0; skip=0

# ── helpers ──────────────────────────────────────────────────────────────────

label_fmt="%-42s"

pass_line()  { printf "PASS  $label_fmt\n" "$1"; }
fail_line()  { printf "FAIL  $label_fmt  %s\n" "$1" "$2"; }
skip_line()  { printf "SKIP  $label_fmt  %s\n" "$1" "$2"; }

indent() { sed 's/^/        /'; }

run_with_timeout() {
    local timeout="$1" outfile="$2"
    shift 2

    "$@" >"$outfile" 2>&1 &
    local pid=$!
    local elapsed=0

    while kill -0 "$pid" 2>/dev/null; do
        if [ "$elapsed" -ge "$timeout" ]; then
            kill "$pid" 2>/dev/null || true
            wait "$pid" 2>/dev/null || true
            return 124
        fi
        sleep 1
        elapsed=$((elapsed+1))
    done

    wait "$pid"
}

# Compare two OUTMC files produced by mcsod.
#   $1 = reference  $2 = generated  $3 = energy tolerance (eV)
# An exact byte match passes immediately (strongest check, and what the
# bit-identical delta path produces on the build machine). Otherwise fall back
# to a tolerant comparison of the energy results (E_min / E_ave / E_std): the
# delta-energy MC path sums the cluster expansion in a different order from the
# full recompute, which on a different compiler/arch can nudge a near-threshold
# Metropolis decision and perturb the trajectory. The sampled energetics must
# still agree to within tolerance. Returns 0 on match, 1 otherwise.
mc_outmc_match() {
    local ref="$1" gen="$2" tol="$3"
    diff -q "$ref" "$gen" >/dev/null 2>&1 && return 0
    local key val_ref val_gen d ok=1
    for key in E_min E_ave E_std; do
        val_ref=$(awk -v k="$key" 'index($0,k){for(i=1;i<=NF;i++) if($i=="="){print $(i+1); exit}}' "$ref")
        val_gen=$(awk -v k="$key" 'index($0,k){for(i=1;i<=NF;i++) if($i=="="){print $(i+1); exit}}' "$gen")
        if [ -z "$val_ref" ] || [ -z "$val_gen" ]; then
            echo "        could not parse $key from OUTMC"; ok=0; break
        fi
        d=$(awk -v a="$val_ref" -v b="$val_gen" 'BEGIN{d=a-b; if(d<0)d=-d; print d}')
        if ! awk -v d="$d" -v t="$tol" 'BEGIN{exit !(d<=t)}'; then
            echo "        $key: |$val_ref - $val_gen| = $d eV > $tol eV"; ok=0
        fi
    done
    [ $ok -eq 1 ] && return 0 || return 1
}

# Compare two mcstatsod thermodynamics.dat files.
#   $1 = reference  $2 = generated
# Exact match passes immediately. Otherwise compare numerically: comment/header
# lines exactly, and each data row's E/F within 0.05 eV and S within 5e-4 eV/K.
# Like the OUTMC check, this tolerates benign MC-trajectory differences (the
# delta-energy + complement-set swap move changes the RNG consumption, so the
# fixed-seed walk explores a different — but statistically equivalent — path).
# Returns 0 on match, 1 otherwise.
thermo_match() {
    local ref="$1" gen="$2"
    diff -q "$ref" "$gen" >/dev/null 2>&1 && return 0
    awk '
      function abs(x){ return x<0 ? -x : x }
      NR==FNR { r[FNR]=$0; nref=FNR; next }
      {
        ln=FNR
        if ($0 ~ /^#/ || $0 ~ /T\/K/) {        # header / comment: exact
          if (r[ln] != $0) { print "        header line "ln" differs"; bad=1 }
          next
        }
        n=split(r[ln], a); split($0, b)
        if (n==0) next
        if (a[1] != b[1]) { print "        T label "ln": "a[1]" vs "b[1]; bad=1; next }
        if (abs(a[2]-b[2]) > 0.05) { print "        E "ln": "a[2]" vs "b[2]; bad=1 }
        if (a[3]=="-" || b[3]=="-") { if (a[3]!=b[3]) { print "        F "ln; bad=1 } }
        else if (abs(a[3]-b[3]) > 0.05) { print "        F "ln": "a[3]" vs "b[3]; bad=1 }
        if (abs(a[4]-b[4]) > 5e-4) { print "        S "ln": "a[4]" vs "b[4]; bad=1 }
      }
      END {
        if (FNR < nref) { print "        generated file shorter than reference ("FNR" vs "nref" lines)"; bad=1 }
        exit bad ? 1 : 0
      }
    ' "$ref" "$gen"
}

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

# ── test_genersod_multitarget ─────────────────────────────────────────────────
# Structure generation with three or more target sites.
#
# Until 0.90 every writer resolved a target-site atom as "sptarget(1), else it
# must be sptarget(2)", so with three targets the third one's atoms were written
# with the second target's symbols — silently, and with no committed example to
# expose it, since example13 and example14 both use FILER = -1.  This test drives
# example13 (Sr on La, Mn on Fe, vacancy on O) at FILER = 0 and diffs the first
# configuration against a committed reference.
#
# $1 = display label   $2 = example directory (3-target INSOD, FILER may be -1)
# $3 = generated structure path relative to the temp project
# $4 = committed reference file
test_genersod_multitarget() {
    local label="$1" dir="$2" struct_rel="$3" ref="$4"

    if [ ! -f "$dir/INSOD" ] || [ ! -f "$dir/SGO" ]; then
        skip_line "$label" "(missing INSOD or SGO)"
        skip=$((skip+1)); return
    fi
    if [ ! -f "$ref" ]; then
        skip_line "$label" "(no committed reference)"
        skip=$((skip+1)); return
    fi

    local tmp; tmp=$(mktemp -d)
    cp "$dir/INSOD" "$dir/SGO" "$tmp/"
    # Ask for CIF output; the example itself ships FILER = -1 (enumeration only).
    sed -i.bak 's/^-1$/0/' "$tmp/INSOD" && rm -f "$tmp/INSOD.bak"

    local out; out=$(cd "$tmp" && PATH="$BIN:$PATH" combsod 2>&1)
    if [ $? -ne 0 ]; then
        fail_line "$label" "[combsod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    out=$(cd "$tmp" && PATH="$BIN:$PATH" genersod 2>&1)
    if [ $? -ne 0 ]; then
        fail_line "$label" "[genersod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    if [ ! -f "$tmp/$struct_rel" ]; then
        fail_line "$label" "[$struct_rel not generated]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    if ! diff -q "$ref" "$tmp/$struct_rel" >/dev/null 2>&1; then
        fail_line "$label" "[$struct_rel differs from reference]"
        diff "$ref" "$tmp/$struct_rel" | head -8 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    rm -rf "$tmp"
    pass_line "$label"; pass=$((pass+1))
}

# ── test_genersod_molecule_second_target ──────────────────────────────────────
# A molecule (@NAME) substituted on the *second* target site, written by LAMMPS.
#
# The LAMMPS writer used to test only the first target's molecule flags, so a
# molecule on a second target was emitted as its raw symbol truncated to three
# characters instead of being expanded into atoms.  Derived from example08 by
# adding a second target (I on Pb) ahead of the existing @MA on Cs, so no new
# type mapping is needed.  Molecule orientations are random, so the check is on
# composition, which is orientation-independent:
#
#   64 Cs - 2 replaced by MA        = 62 Cs   (type 1)
#   64 Pb - 1 replaced by I         = 63 Pb   (type 2)
#   192 I + 1 on the Pb site        = 193 I   (type 3)
#   2 MA = CH3NH3                   = 2 C, 2 N, 12 H  (types 4, 5, 6)
#                                    334 atoms in total
#
# $1 = display label   $2 = example08 directory
test_genersod_molecule_second_target() {
    local label="$1" dir="$2"

    if [ ! -f "$dir/INSOD" ] || [ ! -f "$dir/MA.xyz" ] || [ ! -f "$dir/template_in.lammps" ]; then
        skip_line "$label" "(missing example08 inputs)"
        skip=$((skip+1)); return
    fi

    local tmp; tmp=$(mktemp -d)
    cp "$dir/INSOD" "$dir/SGO" "$dir/MA.xyz" "$dir/template_in.lammps" "$tmp/"

    # sptarget 1 -> "2 1" (Pb then Cs); nsubs 2 -> "1" then "2";
    # newsymbol "@MA Cs" -> "I Pb" then "@MA Cs".
    python3 - "$tmp/INSOD" <<'EOF'
import io, re, sys
p = sys.argv[1]
s = io.open(p).read()
s = re.sub(r"(# sptarget[^\n]*\n)1\n",              r"\g<1>2 1\n",        s)
s = re.sub(r"(# nsubs[^\n]*\n(?:#[^\n]*\n)*)2\n", r"\g<1>1\n2\n",     s)
s = re.sub(r"(# newsymbol[^\n]*\n(?:#[^\n]*\n)*)@MA Cs\n",
           r"\g<1>I Pb\n@MA Cs\n", s)
io.open(p, "w").write(s)
EOF

    local out; out=$(cd "$tmp" && PATH="$BIN:$PATH" combsod 2>&1)
    if [ $? -ne 0 ]; then
        fail_line "$label" "[combsod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    out=$(cd "$tmp" && PATH="$BIN:$PATH" genersod 2>&1)
    if [ $? -ne 0 ]; then
        fail_line "$label" "[genersod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local conf; conf=$(ls "$tmp"/n*/c*/conf.data 2>/dev/null | head -1)
    if [ -z "$conf" ]; then
        fail_line "$label" "[no conf.data generated]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local natoms; natoms=$(awk '/atoms$/ {print $1; exit}' "$conf")
    if [ "$natoms" != "334" ]; then
        fail_line "$label" "[expected 334 atoms, got ${natoms:-none}]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local counts; counts=$(awk '/^Atoms/{f=1;next} f&&NF>=6{print $3}' "$conf" \
                           | sort -n | uniq -c | awk '{printf "%s:%s ", $2, $1}')
    if [ "$counts" != "1:62 2:63 3:193 4:2 5:2 6:12 " ]; then
        fail_line "$label" "[type counts wrong: $counts]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    rm -rf "$tmp"
    pass_line "$label"; pass=$((pass+1))
}

# ── test_gener_best ───────────────────────────────────────────────────────────
# sod_gener.sh -choose <label> resolves a rule to a configuration index.
#
# The workflow it replaces was "inspect OUTSQS and pick the best configuration,
# then sod_gener.sh -choose <index>".  OUTSQS opens with a Rank column followed
# by a Config column, and rank 1 is almost never configuration 1, so the manual
# step invites -choose 1: a valid structure that is not the SQS, with nothing to
# signal the mistake.  The test pins the behaviour that prevents it -- that
# bestSQS resolves to the *second* column -- by requiring its output to equal
# -choose <that column> exactly, and by requiring the two to differ from
# -choose 1 (which would otherwise be a vacuous comparison).
#
# $1 = display label   $2 = example16 directory   $3 = level name
test_gener_best() {
    local label="$1" dir="$2" level="$3"

    if [ ! -f "$dir/INSOD" ] || [ ! -f "$dir/$level/ENSEMBLE" ]; then
        skip_line "$label" "(missing example16 inputs)"
        skip=$((skip+1)); return
    fi

    local tmp; tmp=$(mktemp -d)
    cp "$dir/INSOD" "$dir/SGO" "$dir/EQMATRIX" "$dir/supercell.cif" "$tmp/" 2>/dev/null
    cp "$dir/template_input.gin" "$dir/catlow.lib" "$tmp/" 2>/dev/null
    mkdir -p "$tmp/$level"
    cp "$dir/$level/ENSEMBLE" "$tmp/$level/"
    cat > "$tmp/INSQS" <<'INSQSEOF'
# Maximum cluster order
2

# Cutoff radii (Angstroms) for orders 2..MaxOrder
8.0

# Weights for orders 2..MaxOrder
1.0

# omega and eps_tol for van de Walle scoring
10  1.0E-6
INSQSEOF

    local out rc
    out=$(cd "$tmp/$level" && PATH="$BIN:$PATH" "$BIN/sod_sqs.sh" 2>&1); rc=$?
    if [ $rc -ne 0 ] || [ ! -f "$tmp/$level/OUTSQS" ]; then
        fail_line "$label" "[sod_sqs.sh did not produce OUTSQS]"
        echo "$out" | tail -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local cfg; cfg=$(grep -v '^#' "$tmp/$level/OUTSQS" | head -1 | awk '{print $2}')
    # Total configurations, from the "# Top N of M configurations" header line.
    local nic_count hi_idx
    nic_count=$(awk '/^# Top /{for(i=1;i<=NF;i++) if($i=="of"){print $(i+1)+0; exit}}' \
                "$tmp/$level/OUTSQS")
    if [ -z "$nic_count" ] || [ "$nic_count" -lt 20 ]; then
        skip_line "$label" "(could not read configuration count from OUTSQS)"
        skip=$((skip+1)); rm -rf "$tmp"; return
    fi
    hi_idx=$(( nic_count - 1 ))
    if [ -z "$cfg" ] || [ "$cfg" = "1" ]; then
        skip_line "$label" "(rank 1 is configuration ${cfg:-none}; comparison would be vacuous)"
        skip=$((skip+1)); rm -rf "$tmp"; return
    fi

    out=$(cd "$tmp/$level" && PATH="$BIN:$PATH" "$BIN/sod_gener.sh" -choose bestSQS 2>&1); rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[sod_gener.sh -choose bestSQS exited $rc]"
        echo "$out" | tail -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    if ! echo "$out" | grep -q "selects configuration $cfg"; then
        fail_line "$label" "[bestSQS did not report configuration $cfg]"
        echo "$out" | grep -i best | head -2 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local bestdir; bestdir=$(ls -d "$tmp/$level"/c* 2>/dev/null | head -1)
    if [ -z "$bestdir" ]; then
        fail_line "$label" "[bestSQS generated no configuration directory]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    local keep="$tmp/best_keep"; cp -r "$bestdir" "$keep"
    rm -rf "$tmp/$level"/c*

    # bestSQS must equal -choose <Config>, ...
    out=$(cd "$tmp/$level" && PATH="$BIN:$PATH" "$BIN/sod_gener.sh" -choose "$cfg" 2>&1)
    local choosedir; choosedir=$(ls -d "$tmp/$level"/c* 2>/dev/null | head -1)
    if [ -z "$choosedir" ] || ! diff -r "$keep" "$choosedir" >/dev/null 2>&1; then
        fail_line "$label" "[bestSQS differs from -choose $cfg]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    rm -rf "$tmp/$level"/c*

    # ... and must NOT equal -choose 1, the mistake it exists to prevent.
    out=$(cd "$tmp/$level" && PATH="$BIN:$PATH" "$BIN/sod_gener.sh" -choose 1 2>&1)
    local onedir; onedir=$(ls -d "$tmp/$level"/c* 2>/dev/null | head -1)
    if [ -n "$onedir" ] && diff -r "$keep" "$onedir" >/dev/null 2>&1; then
        fail_line "$label" "[bestSQS matched -choose 1; rank and index are being confused]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Energy labels: plant a known minimum and maximum and require them back.
    rm -rf "$tmp/$level"/c*
    awk -v n="$nic_count" -v hi="$hi_idx" \
        'BEGIN{for(i=1;i<=n;i++){e=-100+(i%97)/97;
                if(i==12)e=-200; if(i==13)e=-199;
                if(i==hi)e=50; if(i==hi-1)e=49;
                printf "%d  %.6f\n", i, e}}' \
        > "$tmp/$level/ENERGIES"

    local lbl want
    for lbl in lowestENERGY:12 highestENERGY:"$hi_idx"; do
        want="${lbl#*:}"; lbl="${lbl%%:*}"
        rm -rf "$tmp/$level"/c*
        out=$(cd "$tmp/$level" && PATH="$BIN:$PATH" "$BIN/sod_gener.sh" -choose "$lbl" 2>&1); rc=$?
        if [ $rc -ne 0 ] || ! echo "$out" | grep -q "selects configuration $want"; then
            fail_line "$label" "[-choose $lbl did not select configuration $want]"
            echo "$out" | grep -i choose | head -2 | indent
            fail=$((fail+1)); rm -rf "$tmp"; return
        fi
    done

    # ... and with a count, the two extremes in order.
    for lbl in "lowestENERGY:12 13" "highestENERGY:$hi_idx $((hi_idx-1))"; do
        want="${lbl#*:}"; lbl="${lbl%%:*}"
        rm -rf "$tmp/$level"/c*
        out=$(cd "$tmp/$level" && PATH="$BIN:$PATH" "$BIN/sod_gener.sh" -choose "$lbl" 2 2>&1); rc=$?
        if [ $rc -ne 0 ] || ! echo "$out" | grep -q "$want"; then
            fail_line "$label" "[-choose $lbl 2 did not select '$want']"
            echo "$out" | grep -iA1 "selects" | head -3 | indent
            fail=$((fail+1)); rm -rf "$tmp"; return
        fi
    done

    # A missing ENERGIES must be an error, and an unknown label must be rejected.
    rm -f "$tmp/$level/ENERGIES"
    out=$(cd "$tmp/$level" && PATH="$BIN:$PATH" "$BIN/sod_gener.sh" -choose lowestENERGY 2>&1); rc=$?
    if [ $rc -eq 0 ]; then
        fail_line "$label" "[-choose lowestENERGY succeeded with no ENERGIES]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    out=$(cd "$tmp/$level" && PATH="$BIN:$PATH" "$BIN/sod_gener.sh" -choose bogus 2>&1); rc=$?
    if [ $rc -eq 0 ] || ! echo "$out" | grep -q "neither a configuration index nor a known selection label"; then
        fail_line "$label" "[unknown label was not rejected]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    # Explicit indices must still work.
    rm -rf "$tmp/$level"/c*
    out=$(cd "$tmp/$level" && PATH="$BIN:$PATH" "$BIN/sod_gener.sh" -choose 3 7 2>&1); rc=$?
    if [ $rc -ne 0 ] || [ "$(ls -d "$tmp/$level"/c* 2>/dev/null | wc -l)" -ne 2 ]; then
        fail_line "$label" "[-choose with two indices no longer generates two configurations]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # A count after the label takes that many, in rank order.  The first three
    # OUTSQS data rows are ranks 1-3, so their Config column is the expectation.
    local want3 got3
    want3=$(grep -v '^#' "$tmp/$level/OUTSQS" | head -3 | awk '{printf "%s ", $2}')
    rm -rf "$tmp/$level"/c*
    out=$(cd "$tmp/$level" && PATH="$BIN:$PATH" "$BIN/sod_gener.sh" -choose bestSQS 3 2>&1); rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[-choose bestSQS 3 exited $rc]"
        echo "$out" | tail -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    got3=$(ls -d "$tmp/$level"/c* 2>/dev/null | sed 's|.*/c||' | sed 's/^0*//' | sort -n | tr '\n' ' ')
    if [ "$got3" != "$(printf '%s ' $(printf '%s\n' $want3 | sort -n))" ]; then
        fail_line "$label" "[-choose bestSQS 3 gave '$got3', expected '$want3']"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Asking for more than OUTSQS lists must fail, and say why.
    out=$(cd "$tmp/$level" && PATH="$BIN:$PATH" "$BIN/sod_gener.sh" -choose bestSQS 100000 2>&1); rc=$?
    if [ $rc -eq 0 ] || ! echo "$out" | grep -q "n_top_sqs"; then
        fail_line "$label" "[over-large bestSQS count not rejected with an n_top_sqs hint]"
        echo "$out" | tail -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # A non-positive count is rejected.
    out=$(cd "$tmp/$level" && PATH="$BIN:$PATH" "$BIN/sod_gener.sh" -choose bestSQS 0 2>&1); rc=$?
    if [ $rc -eq 0 ]; then
        fail_line "$label" "[-choose bestSQS 0 was accepted]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Missing OUTSQS must be an error, not a silent fallback.
    rm -f "$tmp/$level/OUTSQS" && rm -rf "$tmp/$level"/c*
    out=$(cd "$tmp/$level" && PATH="$BIN:$PATH" "$BIN/sod_gener.sh" -choose bestSQS 2>&1); rc=$?
    if [ $rc -eq 0 ]; then
        fail_line "$label" "[bestSQS succeeded with no OUTSQS]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    if ! echo "$out" | grep -q "requires"; then
        fail_line "$label" "[unexpected error for missing OUTSQS]"
        echo "$out" | tail -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    rm -rf "$tmp"
    pass_line "$label"; pass=$((pass+1))
}

# ── test_gener_choose_random ──────────────────────────────────────────────────
# sod_gener.sh -choose random N draws N configurations to generate.
#
# Two properties, one structural and one statistical.
#
# Without replacement is a requirement rather than a nicety: each selection
# becomes a directory cNNN, so a repeated draw would collapse into one and the
# run would quietly produce fewer structures than asked for.  The test pins
# that by asking for the whole ensemble, where any repeat must show up as a
# missing configuration: N = nic has to come back as a permutation of 1..nic.
# It also checks a small draw is distinct, in range, and differs between runs
# (nic = 114 and N = 6 makes a coincidental match a ~1e-8 event), and that
# asking for more than the ensemble holds is an error.
#
# The draw is weighted by degeneracy, so that a subset is a fair sample of
# configuration space rather than of the ENSEMBLE's rows.  That is checked on
# the mean degeneracy of the configuration drawn, whose two models are computed
# from the ENSEMBLE itself: sum(d)/nic if rows were equiprobable, sum(d^2)/sum(d)
# under correct weighting.  For example20 those are 271.2 and 320.2.  A single
# draw is far too noisy to separate them -- measured over 40 trials of N=57, the
# top-degeneracy count ranged 27 to 42 against an equiprobable expectation of 28
# -- so the check aggregates 200 single draws and requires the observed mean to
# land above the midpoint, which puts both false-fail and false-pass below 1e-3.
# The draws run with FILER = -1 so no structures are written: 11 ms each.
#
# $1 = display label   $2 = example20 directory
test_gener_choose_random() {
    local label="$1" dir="$2"

    if [ ! -f "$dir/INSOD" ] || [ ! -f "$dir/SGO" ]; then
        skip_line "$label" "(missing example20 inputs)"
        skip=$((skip+1)); return
    fi

    local tmp; tmp=$(mktemp -d)
    cp "$dir/INSOD" "$dir/SGO" "$dir/MA.xyz" "$dir/FA.xyz" "$tmp/" 2>/dev/null

    local out rc
    out=$(cd "$tmp" && PATH="$BIN:$PATH" combsod 2>&1); rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[combsod error]"
        echo "$out" | tail -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    local ndir; ndir=$(ls -d "$tmp"/n*/ 2>/dev/null | head -1)
    local nic; nic=$(awk 'NR==1{for(i=1;i<=NF;i++) if($i=="configurations;"){print $(i-1); exit}}' "$ndir/ENSEMBLE")
    if [ -z "$nic" ] || [ "$nic" -lt 20 ]; then
        skip_line "$label" "(could not read configuration count from ENSEMBLE)"
        skip=$((skip+1)); rm -rf "$tmp"; return
    fi

    # Helper: run a draw and echo the generated indices, ascending.
    draw() {
        rm -rf "$ndir"/c*
        (cd "$ndir" && PATH="$BIN:$PATH" "$BIN/sod_gener.sh" -choose random "$1" >/dev/null 2>&1) || return 1
        ls -d "$ndir"/c* 2>/dev/null | sed 's|.*/c||' | sed 's/^0*//' | sort -n | tr '\n' ' '
    }

    local a b
    a=$(draw 6) || { fail_line "$label" "[-choose random 6 failed]"; fail=$((fail+1)); rm -rf "$tmp"; return; }
    if [ "$(echo $a | wc -w)" -ne 6 ]; then
        fail_line "$label" "[-choose random 6 generated $(echo $a | wc -w) configurations]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    if [ "$(echo $a | tr ' ' '\n' | sort -u | wc -l)" -ne 6 ]; then
        fail_line "$label" "[-choose random 6 repeated a configuration: $a]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    if echo "$a" | tr ' ' '\n' | awk -v n="$nic" 'NF && ($1<1 || $1>n){exit 1}'; then :; else
        fail_line "$label" "[-choose random 6 produced an out-of-range index: $a]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    b=$(draw 6) || { fail_line "$label" "[second -choose random 6 failed]"; fail=$((fail+1)); rm -rf "$tmp"; return; }
    if [ "$a" = "$b" ]; then
        fail_line "$label" "[two draws of 6 from $nic were identical; not random]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # N = nic must be a permutation: every configuration exactly once.
    local all expect
    all=$(draw "$nic") || { fail_line "$label" "[-choose random $nic failed]"; fail=$((fail+1)); rm -rf "$tmp"; return; }
    expect=$(seq 1 "$nic" | tr '\n' ' ')
    if [ "$all" != "$expect" ]; then
        fail_line "$label" "[-choose random $nic is not a permutation of 1..$nic]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # More than the ensemble holds is an error.
    rm -rf "$ndir"/c*
    out=$(cd "$ndir" && PATH="$BIN:$PATH" "$BIN/sod_gener.sh" -choose random $((nic + 1)) 2>&1); rc=$?
    if [ $rc -eq 0 ]; then
        fail_line "$label" "[-choose random $((nic+1)) was accepted for an ensemble of $nic]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # The draw must be weighted by degeneracy.  Expectations from the ENSEMBLE.
    local means lo hi mid
    means=$(awk 'NF>2 && $1 ~ /^[0-9]+$/ {n++; s+=$2; s2+=$2*$2}
                 END{if(n>0 && s>0) printf "%.4f %.4f", s/n, s2/s}' "$ndir/ENSEMBLE")
    lo=${means%% *}; hi=${means##* }
    if [ -z "$lo" ] || [ -z "$hi" ] || \
       ! awk -v a="$hi" -v b="$lo" 'BEGIN{exit !(a > b*1.05)}'; then
        skip_line "$label" "(degeneracies too flat to test weighting)"
        skip=$((skip+1)); rm -rf "$tmp"; return
    fi
    mid=$(awk -v a="$lo" -v b="$hi" 'BEGIN{printf "%.4f", (a+b)/2}')

    # FILER = -1: resolve and report the draw without writing any structure.
    sed -i.bak 's/^11$/-1/' "$tmp/INSOD" && rm -f "$tmp/INSOD.bak"
    rm -rf "$ndir"/c*
    local ndraw=200 picks
    picks=$(for _ in $(seq 1 $ndraw); do
                (cd "$ndir" && PATH="$BIN:$PATH" "$BIN/sod_gener.sh" -choose random 1 2>/dev/null) \
                  | awk '/proportional to degeneracy/{getline; print $1}'
            done)
    local got; got=$(echo "$picks" | wc -w)
    if [ "$got" -ne "$ndraw" ]; then
        fail_line "$label" "[expected $ndraw draws, captured $got]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    local obs
    obs=$(awk 'NR==FNR{ if(NF>2 && $1 ~ /^[0-9]+$/) d[$1]=$2; next }
               NF{s+=d[$1]; n++} END{if(n>0) printf "%.4f", s/n}' "$ndir/ENSEMBLE" <(echo "$picks"))
    if ! awk -v o="$obs" -v m="$mid" 'BEGIN{exit !(o > m)}'; then
        fail_line "$label" "[draw looks unweighted: mean degeneracy $obs, expected near $hi, equiprobable would give $lo]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    rm -rf "$tmp"
    pass_line "$label"; pass=$((pass+1))
}

# ── test_sqs_no_ensemble_scored ───────────────────────────────────────────────
# sod_sqs.sh run from SODPROJECT/ when the ensemble is one level down.
#
# sod_random.sh writes to nXX/random/, and the sweep over nXX/ then finds no
# ENSEMBLE, printed one "Skipping" line and exited 0 -- a silent no-op reporting
# success, the same shape as the symlink fault fixed in e4ec9295.  It must now
# fail and say where the ensembles actually are.
#
# The two ways of scoring nothing must be told apart: an ENSEMBLE that was never
# found, and one that was found but that sqssod failed on.  Reporting the first
# when it was the second sends the user after the wrong problem, so both
# messages are checked.
#
# $1 = display label   $2 = example16 directory   $3 = level name
test_sqs_no_ensemble_scored() {
    local label="$1" dir="$2" level="$3"

    if [ ! -f "$dir/$level/ENSEMBLE" ]; then
        skip_line "$label" "(missing example16 ensemble)"
        skip=$((skip+1)); return
    fi

    local tmp; tmp=$(mktemp -d)
    cp "$dir/INSOD" "$dir/SGO" "$dir/EQMATRIX" "$dir/supercell.cif" "$tmp/" 2>/dev/null
    # Put the ensemble one level down, as sod_random.sh would.
    mkdir -p "$tmp/$level/random"
    cp "$dir/$level/ENSEMBLE" "$tmp/$level/random/"
    cat > "$tmp/INSQS" <<'INSQSEOF'
# Maximum cluster order
2

# Cutoff radii (Angstroms) for orders 2..MaxOrder
8.0

# Weights for orders 2..MaxOrder
1.0

# omega and eps_tol for van de Walle scoring
10  1.0E-6
INSQSEOF

    local out rc
    out=$(cd "$tmp" && PATH="$BIN:$PATH" "$BIN/sod_sqs.sh" 2>&1); rc=$?
    if [ $rc -eq 0 ]; then
        fail_line "$label" "[scored nothing but exited 0]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    if ! echo "$out" | grep -q "no ENSEMBLE was found to score"; then
        fail_line "$label" "[unexpected error output]"
        echo "$out" | tail -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    if ! echo "$out" | grep -q "$level/random"; then
        fail_line "$label" "[error did not point at $level/random]"
        echo "$out" | tail -4 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Run from the directory that does hold the ensemble: must succeed.
    out=$(cd "$tmp/$level/random" && PATH="$BIN:$PATH" "$BIN/sod_sqs.sh" 2>&1); rc=$?
    if [ $rc -ne 0 ] || [ ! -f "$tmp/$level/random/OUTSQS" ]; then
        fail_line "$label" "[scoring from $level/random/ failed]"
        echo "$out" | tail -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # An ENSEMBLE that is found but cannot be scored must say so, and must not
    # claim the ensemble is missing: without INSQS anywhere, sqssod fails.
    local tmp2; tmp2=$(mktemp -d)
    cp "$dir/INSOD" "$dir/SGO" "$dir/EQMATRIX" "$dir/supercell.cif" "$tmp2/" 2>/dev/null
    mkdir -p "$tmp2/$level"; cp "$dir/$level/ENSEMBLE" "$tmp2/$level/"
    out=$(cd "$tmp2" && PATH="$BIN:$PATH" "$BIN/sod_sqs.sh" 2>&1); rc=$?
    if [ $rc -eq 0 ]; then
        fail_line "$label" "[unscorable ensemble exited 0]"
        fail=$((fail+1)); rm -rf "$tmp" "$tmp2"; return
    fi
    if ! echo "$out" | grep -q "sqssod failed for every ensemble it found"; then
        fail_line "$label" "[scoring failure was not distinguished from a missing ensemble]"
        echo "$out" | tail -3 | indent
        fail=$((fail+1)); rm -rf "$tmp" "$tmp2"; return
    fi
    rm -rf "$tmp2"

    rm -rf "$tmp"
    pass_line "$label"; pass=$((pass+1))
}

# ── test_eqmatrix_only ────────────────────────────────────────────────────────
# sod_comb.sh -eqmatrix-only on a supercell too large to enumerate.
#
# randomsod -sym on requires EQMATRIX, but the only producer of one is combsod,
# whose enumeration is infeasible at exactly the sizes where random sampling is
# wanted: example18 is C(192,96) ~ 3e17 configurations and combsod dies in the
# allocator.  It dies *after* writing supercell.cif, EQMATRIX and OPERATORS, so
# the files were obtainable, but only as the by-product of a failed run that
# exits 1 and prints "problem too large for SOD".  -eqmatrix-only makes that an
# intentional, successful operation.
#
# Checks that the run succeeds, that the three artefacts are byte-identical to
# the ones committed for example18 (i.e. the flag changes when combsod stops,
# not what it writes), that no ENSEMBLE is produced, and that randomsod -sym on
# then works against the result.  Also checks that the site-count validation
# still fires, since the whole point of stopping after setup_targets rather
# than before it is that a bad INSOD is still rejected.
#
# $1 = display label   $2 = example18 directory
test_eqmatrix_only() {
    local label="$1" dir="$2"

    if [ ! -f "$dir/INSOD" ] || [ ! -f "$dir/SGO" ] || [ ! -f "$dir/EQMATRIX" ]; then
        skip_line "$label" "(missing example18 inputs)"
        skip=$((skip+1)); return
    fi

    local tmp; tmp=$(mktemp -d)
    cp "$dir/INSOD" "$dir/SGO" "$tmp/"

    local out rc
    out=$(cd "$tmp" && PATH="$BIN:$PATH" "$BIN/sod_comb.sh" -eqmatrix-only 2>&1); rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[sod_comb.sh -eqmatrix-only exited $rc]"
        echo "$out" | tail -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local f
    for f in supercell.cif EQMATRIX OPERATORS; do
        if [ ! -s "$tmp/$f" ]; then
            fail_line "$label" "[$f not written]"
            fail=$((fail+1)); rm -rf "$tmp"; return
        fi
        if ! cmp -s "$tmp/$f" "$dir/$f"; then
            fail_line "$label" "[$f differs from the committed example18 copy]"
            fail=$((fail+1)); rm -rf "$tmp"; return
        fi
    done

    if [ -e "$tmp/ENSEMBLE" ] || ls "$tmp"/n*/ENSEMBLE >/dev/null 2>&1; then
        fail_line "$label" "[an ENSEMBLE was produced; enumeration did not stop]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # The EQMATRIX is usable: symmetry-folded random sampling needs one.
    out=$(cd "$tmp" && PATH="$BIN:$PATH" "$BIN/sod_random.sh" -nconf 50 -sym on -seed 12345 2>&1); rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[randomsod -sym on failed against the generated EQMATRIX]"
        echo "$out" | tail -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    if ! ls "$tmp"/n*/random/ENSEMBLE >/dev/null 2>&1; then
        fail_line "$label" "[randomsod wrote no ENSEMBLE]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Stopping after setup_targets, not before it, keeps the site-count check.
    local bad; bad=$(mktemp -d)
    cp "$dir/SGO" "$bad/"
    sed 's/^96$/9999/' "$dir/INSOD" > "$bad/INSOD"
    out=$(cd "$bad" && PATH="$BIN:$PATH" "$BIN/sod_comb.sh" -eqmatrix-only 2>&1); rc=$?
    if [ $rc -eq 0 ]; then
        fail_line "$label" "[nsubs exceeding available sites was not rejected]"
        fail=$((fail+1)); rm -rf "$tmp" "$bad"; return
    fi
    if ! echo "$out" | grep -q "exceeds number of available sites"; then
        fail_line "$label" "[unexpected error for oversized nsubs]"
        echo "$out" | tail -3 | indent
        fail=$((fail+1)); rm -rf "$tmp" "$bad"; return
    fi

    rm -rf "$tmp" "$bad"
    pass_line "$label"; pass=$((pass+1))
}

# ── test_genersod_multinary_mol_vac ───────────────────────────────────────────
# Molecules (@) and vacancies (%) on *multi-nary* targets (nk >= 2).
#
# Until 0.90 the @ and % prefixes were only honoured on binary targets and any
# occurrence on a target with nk >= 2 aborted genersod at startup; the generic
# site resolver added in 0.91 handles every slot of every target.  example20
# puts both on both targets at once: the A site is @FA / %VA / @MA, so even the
# remaining species is a molecule, and the X site is Br / %VX / I, a vacancy
# between two plain atoms.  Composition of the 2x2x2 cell:
#
#   1 FA  = 1 C, 2 N,  6 H        6 MA = 6 C, 6 N, 36 H
#   1 A-site vacancy, 1 X-site vacancy, both contributing no atoms
#   8 Pb (untouched B site), 1 Br, 22 I
#                                  88 atoms in total
#
# Molecule orientations are sampled at random and the VASP writer lists species
# in order of first appearance, so the check is on the label -> count mapping,
# which is invariant, rather than on the file or the species line as written.
#
# $1 = display label   $2 = example20 directory
test_genersod_multinary_mol_vac() {
    local label="$1" dir="$2"

    if [ ! -f "$dir/INSOD" ] || [ ! -f "$dir/SGO" ] || \
       [ ! -f "$dir/MA.xyz" ] || [ ! -f "$dir/FA.xyz" ]; then
        skip_line "$label" "(missing example20 inputs)"
        skip=$((skip+1)); return
    fi

    local tmp; tmp=$(mktemp -d)
    cp "$dir/INSOD" "$dir/SGO" "$dir/MA.xyz" "$dir/FA.xyz" "$tmp/"

    local out
    out=$(cd "$tmp" && PATH="$BIN:$PATH" combsod 2>&1)
    if [ $? -ne 0 ]; then
        fail_line "$label" "[combsod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    out=$(cd "$tmp" && PATH="$BIN:$PATH" genersod 2>&1)
    if [ $? -ne 0 ]; then
        fail_line "$label" "[genersod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local poscar; poscar=$(ls "$tmp"/n*/c*/POSCAR 2>/dev/null | head -1)
    if [ -z "$poscar" ]; then
        fail_line "$label" "[no POSCAR generated]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Lines 6 and 7 are the species labels and their counts.
    local comp; comp=$(awk 'NR==6{split($0,s)} NR==7{for(i=1;i<=NF;i++) print s[i]":"$i; exit}' \
                       "$poscar" | sort | tr '\n' ' ')
    if [ "$comp" != "Br:1 C:7 H:42 I:22 N:8 Pb:8 " ]; then
        fail_line "$label" "[composition wrong: $comp]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local natoms; natoms=$(awk 'NR==7{s=0; for(i=1;i<=NF;i++) s+=$i; print s}' "$poscar")
    if [ "$natoms" != "88" ]; then
        fail_line "$label" "[expected 88 atoms, got ${natoms:-none}]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # The vacancy tokens must never reach the structure file as atom labels.
    if grep -qE '(^|[[:space:]])(%|@|VA|VX)' "$poscar"; then
        fail_line "$label" "[vacancy or molecule token leaked into POSCAR]"
        grep -nE '(^|[[:space:]])(%|@|VA|VX)' "$poscar" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    rm -rf "$tmp"
    pass_line "$label"; pass=$((pass+1))
}

# ── test_wrapper_filer_tail ───────────────────────────────────────────────────
# Ensure the shell wrappers read FILER from the last INSOD data line rather than
# literally the last file line, so trailing comments/blanks do not break them.
#
# $1 = display label   $2 = mode ("comb" or "gener")   $3 = example directory
# $4 = expected output path relative to the temp project
test_wrapper_filer_tail() {
    local label="$1" mode="$2" dir="$3" expect_rel="$4"

    if [ ! -s "$dir/INSOD" ]; then
        skip_line "$label" "(empty INSOD)"
        skip=$((skip+1)); return
    fi

    local tmp; tmp=$(mktemp -d)
    rsync -a \
        --exclude='n*/'        \
        --exclude='EQMATRIX'   \
        --exclude='filer'      \
        --exclude='OPERATORS'  \
        --exclude='supercell.cif' \
        --exclude='job_sender' \
        "$dir/" "$tmp/" 2>/dev/null

    printf '\n# trailing comment after FILER\n\n' >> "$tmp/INSOD"

    local out rc
    if [ "$mode" = "comb" ]; then
        out=$(cd "$tmp" && PATH="$BIN:$PATH" "$BIN/sod_comb.sh" 2>&1)
        rc=$?
        if [ $rc -ne 0 ]; then
            fail_line "$label" "[sod_comb.sh error]"
            echo "$out" | head -3 | indent
            fail=$((fail+1)); rm -rf "$tmp"; return
        fi
    elif [ "$mode" = "gener" ]; then
        out=$(cd "$tmp" && PATH="$BIN:$PATH" combsod 2>&1)
        rc=$?
        if [ $rc -ne 0 ]; then
            fail_line "$label" "[combsod setup error]"
            echo "$out" | head -3 | indent
            fail=$((fail+1)); rm -rf "$tmp"; return
        fi
        out=$(cd "$tmp" && PATH="$BIN:$PATH" "$BIN/sod_gener.sh" 2>&1)
        rc=$?
        if [ $rc -ne 0 ]; then
            fail_line "$label" "[sod_gener.sh error]"
            echo "$out" | head -3 | indent
            fail=$((fail+1)); rm -rf "$tmp"; return
        fi
    else
        fail_line "$label" "[unknown test mode: $mode]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    if [ ! -f "$tmp/$expect_rel" ]; then
        fail_line "$label" "[missing $expect_rel]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    pass_line "$label"; pass=$((pass+1))
    rm -rf "$tmp"
}

# ── test_wrapper_model_flag_error ─────────────────────────────────────────────
# Ensure wrapper-level -model parsing fails fast and returns nonzero instead of
# hanging (sod_mc.sh) or returning success after an argument error (sod_cpme.sh).
#
# $1 = display label   $2 = script path   $3 = example directory
test_wrapper_model_flag_error() {
    local label="$1" script_path="$2" dir="$3"

    local tmp; tmp=$(mktemp -d)
    rsync -a "$dir/" "$tmp/" 2>/dev/null

    local out="$tmp/wrapper.err"
    run_with_timeout 5 "$out" bash -lc "cd \"$tmp\" && PATH=\"$BIN:\$PATH\" \"$script_path\" -model"
    local rc=$?

    if [ $rc -eq 124 ]; then
        fail_line "$label" "[wrapper hung on missing -model argument]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    if [ $rc -eq 0 ]; then
        fail_line "$label" "[wrapper returned success on missing -model argument]"
        sed -n '1,3p' "$out" | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    if ! grep -q "Error: -model requires a filename argument." "$out"; then
        fail_line "$label" "[missing expected error message]"
        sed -n '1,5p' "$out" | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    pass_line "$label"; pass=$((pass+1))
    rm -rf "$tmp"
}

# ── CELL completeness policy ─────────────────────────────────────────────────

test_cell_completeness() {
    local label="$1"
    local tmp source_file config_name out mismatch=0

    # Both halves replay committed calculator output (VASP CONTCAR, GULP
    # output.gout) from example01.  release.sh scrubs those from the public
    # package, so in a released tree there is nothing to replay: skip rather
    # than fail, as for any other missing optional input.
    if ! ls "$EX/example01/FILER11_vasp/n04"/c*/CONTCAR >/dev/null 2>&1 || \
       ! ls "$EX/example01/FILER1_gulp/n04"/c*/output.gout >/dev/null 2>&1; then
        skip_line "$label" "(no committed calculator output — excluded from releases)"
        skip=$((skip+1)); return
    fi

    tmp=$(mktemp -d)

    # Complete VASP results invoked from SODPROJECT/ must write nXX/CELL (not a
    # combined root CELL), preserving the established extractor output exactly.
    cp "$EX/example01/FILER11_vasp/INSOD" "$EX/example01/FILER11_vasp/SGO" "$tmp/"
    mkdir "$tmp/n04"
    cp "$EX/example01/FILER11_vasp/n04/ENSEMBLE" "$tmp/n04/"
    for source_file in "$EX/example01/FILER11_vasp/n04"/c*/CONTCAR; do
        config_name=$(basename "$(dirname "$source_file")")
        mkdir "$tmp/n04/$config_name"
        cp "$source_file" "$tmp/n04/$config_name/"
    done
    out=$(cd "$tmp" && PATH="$BIN:$PATH" sod_vasp_cell.sh cub 2>&1)
    if [ $? -ne 0 ] || [ ! -f "$tmp/n04/CELL" ] || [ -f "$tmp/CELL" ] || \
       ! diff -q "$EX/example01/FILER11_vasp/n04/CELL" "$tmp/n04/CELL" >/dev/null; then
        fail_line "$label" "[complete VASP results did not produce the reference n04/CELL]"
        mismatch=1
    fi

    # Removing one result makes the set sparse: CELL must be removed/skipped
    # rather than shift every later positional row onto the wrong configuration.
    rm -f "$tmp/n04/c71/CONTCAR"
    out=$(cd "$tmp/n04" && PATH="$BIN:$PATH" sod_vasp_cell.sh cub 2>&1)
    if [ $? -ne 0 ] || [ -f "$tmp/n04/CELL" ] || \
       ! printf '%s' "$out" | grep -q "70 of 71 ENSEMBLE configurations"; then
        fail_line "$label" "[sparse VASP results were allowed to write CELL]"
        mismatch=1
    fi

    rm -rf "$tmp"
    tmp=$(mktemp -d)
    cp "$EX/example01/FILER1_gulp/INSOD" "$EX/example01/FILER1_gulp/SGO" "$tmp/"
    mkdir "$tmp/n04"
    cp "$EX/example01/FILER1_gulp/n04/ENSEMBLE" "$tmp/n04/"
    for source_file in "$EX/example01/FILER1_gulp/n04"/c*/output.gout; do
        config_name=$(basename "$(dirname "$source_file")")
        mkdir "$tmp/n04/$config_name"
        cp "$source_file" "$tmp/n04/$config_name/"
    done
    out=$(cd "$tmp/n04" && PATH="$BIN:$PATH" sod_gulp_cell.sh cub 2>&1)
    if [ $? -ne 0 ] || [ ! -f "$tmp/n04/CELL" ] || \
       [ "$(wc -l < "$tmp/n04/CELL")" -ne 71 ]; then
        fail_line "$label" "[complete GULP results did not produce a 71-row CELL]"
        mismatch=1
    fi

    # A spare cYY directory holding no output.gout leaves coverage complete, so
    # the extractor must skip it silently rather than let awk report it fatal.
    mkdir "$tmp/n04/c99"
    out=$(cd "$tmp/n04" && PATH="$BIN:$PATH" sod_gulp_cell.sh cub 2>&1)
    if [ $? -ne 0 ] || [ "$(wc -l < "$tmp/n04/CELL")" -ne 71 ] || \
       printf '%s' "$out" | grep -q "awk"; then
        fail_line "$label" "[a spare cYY directory disturbed GULP CELL extraction]"
        printf '%s\n' "$out" | sed -n '1,3p' | indent
        mismatch=1
    fi
    rmdir "$tmp/n04/c99"

    rm -f "$tmp/n04/c71/output.gout"
    out=$(cd "$tmp/n04" && PATH="$BIN:$PATH" sod_gulp_cell.sh cub 2>&1)
    if [ $? -ne 0 ] || [ -f "$tmp/n04/CELL" ] || \
       ! printf '%s' "$out" | grep -q "70 of 71 ENSEMBLE configurations"; then
        fail_line "$label" "[sparse GULP results were allowed to write CELL]"
        mismatch=1
    fi

    rm -rf "$tmp"
    if [ $mismatch -eq 0 ]; then
        pass_line "$label"; pass=$((pass+1))
    else
        fail=$((fail+1))
    fi
}

# ── MACE result-file protection (preflight only, no MLIP stack needed) ─────

# ── test_mace_private_api ─────────────────────────────────────────────────────
# sod_mace calls a few private nvalchemi-toolkit APIs.  mace_backend.PRIVATE_APIS
# is the authoritative list of them and check_private_api() enforces it before a
# run starts, so a toolkit upgrade that renames or drops one is caught here by
# `make test` rather than by a user, minutes into a relaxation.  The test also
# plants a bogus entry, so a check that silently passed everything would fail.
test_mace_private_api() {
    local label="$1"
    local py="${SOD_PYTHON:-python3}"

    if ! command -v "$py" >/dev/null 2>&1; then
        skip_line "$label" "(no python interpreter — set SOD_PYTHON)"
        skip=$((skip+1)); return
    fi
    if ! "$py" -c 'import torch, ase, mace, nvalchemi' >/dev/null 2>&1; then
        skip_line "$label" "(python MACE stack not installed)"
        skip=$((skip+1)); return
    fi

    local out
    out=$("$py" - "$ROOT/pysod" <<'PYEOF' 2>&1
import sys

sys.path.insert(0, sys.argv[1])
import mace_backend as mb

mb.check_private_api()

# A check that accepted anything would be worse than none: make sure it fails.
mb.PRIVATE_APIS = mb.PRIVATE_APIS + (("FIRE2", mb.FIRE2, "_sod_absent_api"),)
try:
    mb.check_private_api()
except RuntimeError:
    pass
else:
    raise SystemExit("check_private_api accepted a missing attribute")

print("ok")
PYEOF
)
    if [ $? -ne 0 ] || ! printf '%s' "$out" | grep -q "^ok$"; then
        fail_line "$label" "[private API check failed]"
        printf '%s\n' "$out" | head -5 | indent
        fail=$((fail+1)); return
    fi

    pass_line "$label"; pass=$((pass+1))
}

test_mace_result_protection() {
    local label="$1"
    local tmp py out i mismatch=0
    py="${SOD_PYTHON:-python3}"
    if ! command -v "$py" >/dev/null 2>&1; then
        skip_line "$label" "(no python interpreter — set SOD_PYTHON)"
        skip=$((skip+1)); return
    fi

    tmp=$(mktemp -d)
    touch "$tmp/INSOD" "$tmp/SGO"
    mkdir -p "$tmp/n01"
    printf 'Enumerated ensemble: 40 configurations; sum_degeneracies = 40\n' > "$tmp/n01/ENSEMBLE"
    for i in $(seq -w 1 40); do
        mkdir -p "$tmp/n01/c$i"
        touch "$tmp/n01/c$i/relaxed_POSCAR"
    done

    # The clash list must name the total and stay short: a level can hold
    # thousands of cYY/ directories.
    out=$(cd "$tmp" && "$py" "$ROOT/pysod/sod_mace.py" -relax -structure POSCAR 2>&1)
    if [ $? -eq 0 ] || ! printf '%s' "$out" | grep -q "40 existing result files" || \
       ! printf '%s' "$out" | grep -q "and 32 more" || \
       [ "$(printf '%s' "$out" | head -1 | wc -c)" -gt 400 ]; then
        fail_line "$label" "[relaxed-structure clash list was not summarised]"
        printf '%s\n' "$out" | sed -n '1,2p' | cut -c1-120 | indent
        mismatch=1
    fi
    rm -rf "$tmp"

    # A rerun narrower than the one before it must not silently keep the wider
    # run's summary files: they pair with the ENERGIES that produced them.
    tmp=$(mktemp -d)
    touch "$tmp/INSOD" "$tmp/SGO"
    mkdir -p "$tmp/n01/c1"
    printf 'Enumerated ensemble: 1 configurations; sum_degeneracies = 1\n' > "$tmp/n01/ENSEMBLE"
    echo old > "$tmp/n01/MACE_RELAXATION.dat"
    echo old > "$tmp/n01/CELL"
    echo old > "$tmp/n01/ENTHALPIES"

    out=$(cd "$tmp" && "$py" "$ROOT/pysod/sod_mace.py" -structure POSCAR 2>&1)
    if [ $? -eq 0 ]; then
        fail_line "$label" "[single-point rerun ignored stale summary files]"
        mismatch=1
    fi
    for i in ENTHALPIES MACE_RELAXATION.dat CELL; do
        if ! printf '%s' "$out" | grep -q "n01/$i"; then
            fail_line "$label" "[stale $i was not reported by the preflight]"
            mismatch=1
        fi
    done

    # -force removes them before the run proceeds, whether or not the MLIP
    # stack is present for the run itself.
    (cd "$tmp" && "$py" "$ROOT/pysod/sod_mace.py" -structure POSCAR -force >/dev/null 2>&1)
    for i in ENTHALPIES MACE_RELAXATION.dat CELL; do
        if [ -f "$tmp/n01/$i" ]; then
            fail_line "$label" "[-force left stale $i in place]"
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

# ── ENSEMBLE format version 3 is the only accepted format ────────────────────

test_ensemble_v3_only() {
    local label="$1"
    local tmp source_file config_name out rc mismatch=0

    # Replays committed VASP CONTCARs from example01, which release.sh scrubs
    # from the public package; skip in a released tree rather than fail.
    if ! ls "$EX/example01/FILER11_vasp/n04"/c*/CONTCAR >/dev/null 2>&1; then
        skip_line "$label" "(no committed calculator output — excluded from releases)"
        skip=$((skip+1)); return
    fi

    tmp=$(mktemp -d)

    # A complete VASP level, plus a directory whose name starts with c<digit>
    # but is not a configuration directory. It must be ignored, not parsed as a
    # base-10 index, and must not add a spurious CELL row.
    cp "$EX/example01/FILER11_vasp/INSOD" "$EX/example01/FILER11_vasp/SGO" "$tmp/"
    mkdir "$tmp/n04"
    cp "$EX/example01/FILER11_vasp/n04/ENSEMBLE" "$tmp/n04/"
    for source_file in "$EX/example01/FILER11_vasp/n04"/c*/CONTCAR; do
        config_name=$(basename "$(dirname "$source_file")")
        mkdir "$tmp/n04/$config_name"
        cp "$source_file" "$tmp/n04/$config_name/"
    done
    mkdir "$tmp/n04/c01_backup"
    cp "$tmp/n04/c01/CONTCAR" "$tmp/n04/c01_backup/"

    out=$(cd "$tmp/n04" && PATH="$BIN:$PATH" sod_vasp_cell.sh cub 2>&1)
    rc=$?
    if [ $rc -ne 0 ] || [ ! -f "$tmp/n04/CELL" ] || \
       ! diff -q "$EX/example01/FILER11_vasp/n04/CELL" "$tmp/n04/CELL" >/dev/null; then
        fail_line "$label" "[a non-numeric cYY-like directory broke CELL extraction]"
        printf '%s\n' "$out" | sed -n '1,3p' | indent
        mismatch=1
    fi

    # A version 2 ENSEMBLE must be a hard error, and must leave any existing
    # CELL untouched: coverage cannot be judged, so nothing is known to be stale.
    rm -rf "$tmp/n04/c01_backup"
    printf '# SOD ENSEMBLE format version 2\n4 substitutions in 32 sites\n71\n' > "$tmp/n04/ENSEMBLE"
    out=$(cd "$tmp/n04" && PATH="$BIN:$PATH" sod_vasp_cell.sh cub 2>&1)
    rc=$?
    if [ $rc -eq 0 ] || [ ! -f "$tmp/n04/CELL" ] || \
       ! printf '%s' "$out" | grep -q "is not a version 3 ENSEMBLE"; then
        fail_line "$label" "[a version 2 ENSEMBLE did not fail loudly and preserve CELL]"
        printf '%s\n' "$out" | sed -n '1,3p' | indent
        mismatch=1
    fi
    rm -rf "$tmp"

    # statsod reads ENSEMBLE with its own parser; it must reject version 2 too.
    tmp=$(mktemp -d)
    printf '# SOD ENSEMBLE format version 2\n2 substitutions in 6 sites\n3\n1 3 1 2\n2 6 1 3\n3 3 2 5\n' > "$tmp/ENSEMBLE"
    printf '1 -10.0\n2 -11.0\n3 -12.0\n' > "$tmp/ENERGIES"
    out=$(cd "$tmp" && "$BIN/statsod" 2>&1)
    rc=$?
    if [ $rc -eq 0 ] || ! printf '%s' "$out" | grep -q "not a version 3"; then
        fail_line "$label" "[statsod accepted a version 2 ENSEMBLE]"
        printf '%s\n' "$out" | sed -n '1,4p' | indent
        mismatch=1
    fi
    rm -rf "$tmp"

    # invertENSEMBLE writes version 3, so inverting twice returns the original.
    tmp=$(mktemp -d)
    printf 'Enumerated ensemble: 3 configurations; sum_degeneracies = 12\n6 Fe sites -> 2 Al 4 Fe\n# Configuration  Degeneracy  Al_positions\n1 3 1 2\n2 6 1 3\n3 3 2 5\n' > "$tmp/ENSEMBLE_original"
    cp "$tmp/ENSEMBLE_original" "$tmp/reference"
    out=$(cd "$tmp" && "$BIN/invertENSEMBLE" 2>&1)
    if [ $? -ne 0 ] || ! grep -q "^6 Al sites -> 4 Fe 2 Al$" "$tmp/ENSEMBLE_inverted"; then
        fail_line "$label" "[invertENSEMBLE did not write a version 3 target line]"
        printf '%s\n' "$out" | sed -n '1,4p' | indent
        mismatch=1
    fi
    cp "$tmp/ENSEMBLE_inverted" "$tmp/ENSEMBLE_original"
    (cd "$tmp" && "$BIN/invertENSEMBLE" >/dev/null 2>&1)
    if ! diff -q "$tmp/reference" "$tmp/ENSEMBLE_inverted" >/dev/null; then
        fail_line "$label" "[inverting twice did not reproduce the original ENSEMBLE]"
        diff "$tmp/reference" "$tmp/ENSEMBLE_inverted" | sed -n '1,6p' | indent
        mismatch=1
    fi

    printf '# SOD ENSEMBLE format version 2\n2 substitutions in 6 sites\n3\n1 3 1 2\n' > "$tmp/ENSEMBLE_original"
    out=$(cd "$tmp" && "$BIN/invertENSEMBLE" 2>&1)
    if [ $? -eq 0 ] || ! printf '%s' "$out" | grep -q "not a version 3 ENSEMBLE"; then
        fail_line "$label" "[invertENSEMBLE accepted a version 2 ENSEMBLE]"
        mismatch=1
    fi
    rm -rf "$tmp"

    if [ $mismatch -eq 0 ]; then
        pass_line "$label"; pass=$((pass+1))
    else
        fail=$((fail+1))
    fi
}

test_energy_enthalpy_extractors() {
    local label="$1"
    local tmp got expected out mismatch=0
    tmp=$(mktemp -d)
    touch "$tmp/INSOD" "$tmp/SGO"

    mkdir -p "$tmp/n01/c02" "$tmp/n01/c17"
    cat > "$tmp/n01/c02/OUTCAR" <<'EOF'
  energy  without entropy=      -12.2000  energy(sigma->0) =      -12.0000
  enthalpy is  TOTEN     =      -11.8000 eV   P V=       0.2000
  energy  without entropy=      -10.1000  energy(sigma->0) =      -10.0000
  enthalpy is  TOTEN     =       -9.6000 eV   P V=       0.4000
EOF
    cat > "$tmp/n01/c17/OUTCAR" <<'EOF'
  energy  without entropy=      -11.1000  energy(sigma->0) =      -11.0000
EOF

    (cd "$tmp/n01" && PATH="$BIN:$PATH" sod_vasp_ener.sh >/dev/null 2>&1)
    got=$(cat "$tmp/n01/ENERGIES" 2>/dev/null)
    expected=$(printf "2  -10.0000\n17  -11.0000")
    if [ "$got" != "$expected" ]; then
        fail_line "$label" "[VASP ENERGIES did not keep final energy(sigma->0)]"
        mismatch=1
    fi

    # c17 has no "P V=" line, so it carries no enthalpy: it is missing from
    # ENTHALPIES rather than recorded with its internal energy.
    (cd "$tmp/n01" && PATH="$BIN:$PATH" sod_vasp_enth.sh >/dev/null 2>&1)
    got=$(cat "$tmp/n01/ENTHALPIES" 2>/dev/null)
    expected="2  -9.6000000000"
    if [ "$got" != "$expected" ]; then
        fail_line "$label" "[VASP ENTHALPIES did not compute E+PV with sparse indices]"
        mismatch=1
    fi

    # A level with no PSTRESS anywhere must write no ENTHALPIES and clear a
    # stale one, matching sod_mace.sh at zero pressure.
    mkdir -p "$tmp/n03/c01"
    cp "$tmp/n01/c17/OUTCAR" "$tmp/n03/c01/OUTCAR"
    echo "stale" > "$tmp/n03/ENTHALPIES"
    out=$(cd "$tmp/n03" && PATH="$BIN:$PATH" sod_vasp_enth.sh 2>&1)
    if [ -f "$tmp/n03/ENTHALPIES" ] || ! printf '%s' "$out" | grep -q "zero pressure"; then
        fail_line "$label" "[zero-pressure VASP level still produced ENTHALPIES]"
        mismatch=1
    fi

    mkdir -p "$tmp/n02/c02" "$tmp/n02/c17"
    cat > "$tmp/n02/c02/output.gout" <<'EOF'
  Final enthalpy =        -24.000000 eV
  Pressure*volume =         1.000000 eV
  Final enthalpy =        -20.000000 eV
  Pressure*volume =         1.500000 eV
EOF
    cat > "$tmp/n02/c17/output.gout" <<'EOF'
  Final energy =        -30.000000 eV
EOF

    (cd "$tmp/n02" && PATH="$BIN:$PATH" sod_gulp_ener.sh >/dev/null 2>&1)
    got=$(cat "$tmp/n02/ENERGIES" 2>/dev/null)
    expected=$(printf "2  -21.5000000000\n17  -30.0000000000")
    if [ "$got" != "$expected" ]; then
        fail_line "$label" "[GULP ENERGIES did not subtract Pressure*volume from enthalpy]"
        mismatch=1
    fi

    # As for VASP, c17 has no Pressure*volume term and so has no enthalpy.
    (cd "$tmp/n02" && PATH="$BIN:$PATH" sod_gulp_enth.sh >/dev/null 2>&1)
    got=$(cat "$tmp/n02/ENTHALPIES" 2>/dev/null)
    expected="2  -20.0000000000"
    if [ "$got" != "$expected" ]; then
        fail_line "$label" "[GULP ENTHALPIES did not keep native final enthalpy]"
        mismatch=1
    fi

    mkdir -p "$tmp/n04/c01"
    cp "$tmp/n02/c17/output.gout" "$tmp/n04/c01/output.gout"
    echo "stale" > "$tmp/n04/ENTHALPIES"
    out=$(cd "$tmp/n04" && PATH="$BIN:$PATH" sod_gulp_enth.sh 2>&1)
    if [ -f "$tmp/n04/ENTHALPIES" ] || ! printf '%s' "$out" | grep -q "zero pressure"; then
        fail_line "$label" "[zero-pressure GULP level still produced ENTHALPIES]"
        mismatch=1
    fi

    rm -rf "$tmp"
    if [ $mismatch -eq 0 ]; then
        pass_line "$label"; pass=$((pass+1))
    else
        fail=$((fail+1))
    fi
}

# ── test_stat ────────────────────────────────────────────────────────────────
# $1 = display label   $2 = n-level directory (must contain ENSEMBLE, ENERGIES;
#      optionally TEMPERATURES, DATA and SPECTRA)
# If TEMPERATURES is present it is used; otherwise statsod's default 0/300/1000 K.
# Reference files that must already be committed in $2:
#   thermodynamics.dat   (always)
#   ave_data.dat         (if DATA is present)
#   ave_spectra.dat      (if SPECTRA is present)
test_sod_mace() {
    local label="$1" xdir="$2" nxx="$3"
    local py="${SOD_PYTHON:-python3}"

    if [ ! -f "$xdir/INSOD" ] || [ ! -f "$xdir/SGO" ]; then
        skip_line "$label" "(missing INSOD or SGO)"
        skip=$((skip+1)); return
    fi
    # The MLIP stack is an optional dependency: CI installs gfortran only, so a
    # missing torch/MACE/ALCHEMI must skip rather than fail. Set SOD_PYTHON to
    # the interpreter that has them (see pysod/README.md).
    if ! command -v "$py" >/dev/null 2>&1; then
        skip_line "$label" "(no python interpreter — set SOD_PYTHON)"
        skip=$((skip+1)); return
    fi
    if ! "$py" -c 'import torch, ase, mace, nvalchemi' >/dev/null 2>&1; then
        skip_line "$label" "(python MACE stack not installed)"
        skip=$((skip+1)); return
    fi

    # Only a handful of configurations: this checks identifier mapping and output
    # format, not energies, which depend on model version, device and precision.
    local tmp; tmp=$(mktemp -d)
    cp "$xdir/INSOD" "$xdir/SGO" "$tmp/"
    mkdir -p "$tmp/$nxx"
    local n=0 c
    for c in "$xdir/$nxx"/c*/; do
        [ -d "$c" ] || continue
        cp -r "$c" "$tmp/$nxx/" || true
        # Drop any relaxed structure a previous sod_mace.sh run left in the
        # example tree: copying it in would trip the result-file protection.
        rm -f "$tmp/$nxx/$(basename "${c%/}")"/relaxed.* \
              "$tmp/$nxx/$(basename "${c%/}")"/relaxed_*
        n=$((n+1))
        [ $n -ge 4 ] && break
    done
    if [ $n -eq 0 ]; then
        skip_line "$label" "(no cYY/ folders in $nxx)"
        skip=$((skip+1)); rm -rf "$tmp"; return
    fi

    local out rc
    out=$(cd "$tmp/$nxx" && PATH="$BIN:$PATH" sod_mace.sh -device cpu -cueq off -batch 2 -q 2>&1)
    rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[sod_mace.sh error]"
        echo "$out" | tail -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local mismatch=0
    if [ ! -f "$tmp/$nxx/ENERGIES" ]; then
        fail_line "$label" "[no ENERGIES written]"
        mismatch=1
    else
        # Two-column "m  E" with dense 1-based indices, as read_energies_file expects.
        local check
        check=$(awk -v want="$n" '
            /^[[:space:]]*(#|$)/ { next }
            { rows++
              if ($1 != rows) { printf "index %s out of order at row %d\n", $1, rows; exit }
              if ($2 !~ /^-?[0-9]+\.?[0-9]*([eE][-+]?[0-9]+)?$/) { printf "bad energy %s\n", $2; exit }
              if (NF != 2) { printf "expected 2 columns, got %d\n", NF; exit } }
            END { if (rows != want) printf "expected %d rows, got %d\n", want, rows }
        ' "$tmp/$nxx/ENERGIES")
        if [ -n "$check" ]; then
            fail_line "$label" "[$check]"
            mismatch=1
        fi
    fi

    # A second run must refuse to overwrite existing energies unless -force.
    local before after
    before=$(cat "$tmp/$nxx/ENERGIES" 2>/dev/null)
    if (cd "$tmp/$nxx" && PATH="$BIN:$PATH" sod_mace.sh -device cpu -cueq off -q >/dev/null 2>&1); then
        fail_line "$label" "[rerun overwrote ENERGIES without -force]"
        mismatch=1
    fi
    after=$(cat "$tmp/$nxx/ENERGIES" 2>/dev/null)
    if [ "$before" != "$after" ]; then
        fail_line "$label" "[ENERGIES changed by a refused rerun]"
        mismatch=1
    fi

    # A stale pressure-run ENTHALPIES file must not be left beside fresh
    # zero-pressure energies. Without -force the run refuses; with -force it
    # removes the stale file and proceeds.
    rm -f "$tmp/$nxx/ENERGIES"
    printf "1  0.0\n" > "$tmp/$nxx/ENTHALPIES"
    if (cd "$tmp/$nxx" && PATH="$BIN:$PATH" sod_mace.sh -device cpu -cueq off -q >/dev/null 2>&1); then
        fail_line "$label" "[zero-pressure rerun ignored stale ENTHALPIES]"
        mismatch=1
    fi
    out=$(cd "$tmp/$nxx" && PATH="$BIN:$PATH" sod_mace.sh -device cpu -cueq off -force -q 2>&1)
    rc=$?
    if [ $rc -ne 0 ] || [ -f "$tmp/$nxx/ENTHALPIES" ]; then
        fail_line "$label" "[-force did not remove stale ENTHALPIES]"
        [ -n "$out" ] && echo "$out" | tail -3 | indent
        mismatch=1
    fi

    rm -rf "$tmp"
    if [ $mismatch -eq 0 ]; then
        pass_line "$label"; pass=$((pass+1))
    else
        fail=$((fail+1))
    fi
}

test_sod_mace_relax() {
    local label="$1" xdir="$2" nxx="$3"
    local py="${SOD_PYTHON:-python3}"

    if [ ! -f "$xdir/INSOD" ] || [ ! -f "$xdir/SGO" ]; then
        skip_line "$label" "(missing INSOD or SGO)"
        skip=$((skip+1)); return
    fi
    if ! command -v "$py" >/dev/null 2>&1 || \
       ! "$py" -c 'import torch, ase, mace, nvalchemi' >/dev/null 2>&1; then
        skip_line "$label" "(python MACE stack not installed)"
        skip=$((skip+1)); return
    fi

    local tmp; tmp=$(mktemp -d)
    cp "$xdir/INSOD" "$xdir/SGO" "$tmp/"
    mkdir -p "$tmp/$nxx"
    local n=0 c
    for c in "$xdir/$nxx"/c*/; do
        [ -d "$c" ] || continue
        cp -r "$c" "$tmp/$nxx/" || true
        # Drop any relaxed structure a previous sod_mace.sh run left in the
        # example tree: copying it in would trip the result-file protection.
        rm -f "$tmp/$nxx/$(basename "${c%/}")"/relaxed.* \
              "$tmp/$nxx/$(basename "${c%/}")"/relaxed_*
        n=$((n+1))
        [ $n -ge 4 ] && break
    done

    # batch 2 over 4 structures guarantees the refill path actually fires.
    local mode rc out
    for mode in fixed refill; do
        out=$(cd "$tmp/$nxx" && PATH="$BIN:$PATH" sod_mace.sh -device cpu -cueq off \
              -batch 2 -relax -batchmode "$mode" -maxsteps 20 -force -q 2>&1)
        rc=$?
        if [ $rc -ne 0 ]; then
            fail_line "$label" "[-batchmode $mode error]"
            echo "$out" | tail -3 | indent
            fail=$((fail+1)); rm -rf "$tmp"; return
        fi
        cp "$tmp/$nxx/MACE_RELAXATION.dat" "$tmp/table_$mode.dat"
    done

    local mismatch=0
    # Relaxed geometries are written beside their inputs, never over them.
    local missing=0 i
    for i in $(seq 1 "$n"); do
        local d; d=$(printf "%s/%s/c%02d" "$tmp" "$nxx" "$i")
        [ -f "$d/relaxed.cif" ] || missing=$((missing+1))
        [ -f "$d/configuration.cif" ] || missing=$((missing+1))
    done
    if [ $missing -ne 0 ]; then
        fail_line "$label" "[$missing expected cif file(s) missing]"
        mismatch=1
    fi

    # Refill reuses batch slots, so a broken source_index mapping would attach
    # results to the wrong configuration. Configurations here differ by ~0.1 eV,
    # far more than the ~1e-4 eV float32 spread between the two paths.
    local check
    check=$(awk '
        FNR==NR { if (!/^#/) { e[$1]=$3; s[$1]=$5; n++ } ; next }
        !/^#/ {
            if (!($1 in e)) { printf "config %s only in refill\n", $1; exit }
            d = e[$1] - $3; if (d < 0) d = -d
            if (d > 1.0e-3) { printf "config %s energy differs by %.3e eV\n", $1, d; exit }
            if (s[$1] != $5) { printf "config %s steps %s vs %s\n", $1, s[$1], $5; exit }
            m++
        }
        END { if (m != n) printf "fixed has %d rows, refill %d\n", n, m }
    ' "$tmp/table_fixed.dat" "$tmp/table_refill.dat")
    if [ -n "$check" ]; then
        fail_line "$label" "[$check]"
        mismatch=1
    fi

    # Removing ENERGIES used to let a relaxation rerun overwrite the remaining
    # relaxation artifacts. The preflight must still catch MACE_RELAXATION.dat
    # or cYY/relaxed.cif before doing any work.
    rm -f "$tmp/$nxx/ENERGIES"
    before=$(cat "$tmp/$nxx/MACE_RELAXATION.dat" 2>/dev/null)
    out=$(cd "$tmp/$nxx" && PATH="$BIN:$PATH" sod_mace.sh -device cpu -cueq off \
          -batch 2 -relax -maxsteps 20 -q 2>&1)
    rc=$?
    after=$(cat "$tmp/$nxx/MACE_RELAXATION.dat" 2>/dev/null)
    if [ $rc -eq 0 ] || ! printf '%s' "$out" | grep -q "MACE_RELAXATION.dat"; then
        fail_line "$label" "[rerun did not refuse existing relaxation artifacts]"
        mismatch=1
    fi
    if [ "$before" != "$after" ]; then
        fail_line "$label" "[MACE_RELAXATION.dat changed by a refused rerun]"
        mismatch=1
    fi

    # -writerelaxed no: same energies and step counts, no per-configuration
    # structure written. For runs with too many configurations to keep one file
    # each, so the level summaries must still be complete.
    rm -f "$tmp/$nxx/ENERGIES" "$tmp/$nxx/MACE_RELAXATION.dat"
    for i in $(seq 1 "$n"); do
        rm -f "$(printf "%s/%s/c%02d" "$tmp" "$nxx" "$i")/relaxed.cif"
    done
    out=$(cd "$tmp/$nxx" && PATH="$BIN:$PATH" sod_mace.sh -device cpu -cueq off \
          -batch 2 -relax -maxsteps 20 -writerelaxed no -q 2>&1)
    rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[-writerelaxed no error]"
        echo "$out" | tail -3 | indent
        mismatch=1
    else
        local wrote=0
        for i in $(seq 1 "$n"); do
            [ -f "$(printf "%s/%s/c%02d" "$tmp" "$nxx" "$i")/relaxed.cif" ] && wrote=$((wrote+1))
        done
        if [ $wrote -ne 0 ]; then
            fail_line "$label" "[-writerelaxed no wrote $wrote relaxed structure(s)]"
            mismatch=1
        fi
        if [ ! -f "$tmp/$nxx/ENERGIES" ] || [ ! -f "$tmp/$nxx/MACE_RELAXATION.dat" ]; then
            fail_line "$label" "[-writerelaxed no skipped a level summary file]"
            mismatch=1
        else
            check=$(awk '
                FNR==NR { if (!/^#/) { e[$1]=$3; s[$1]=$5; n++ } ; next }
                !/^#/ {
                    d = e[$1] - $3; if (d < 0) d = -d
                    if (d > 1.0e-3) { printf "config %s energy differs by %.3e eV\n", $1, d; exit }
                    if (s[$1] != $5) { printf "config %s steps %s vs %s\n", $1, s[$1], $5; exit }
                    m++
                }
                END { if (m != n) printf "fixed has %d rows, writerelaxed-no %d\n", n, m }
            ' "$tmp/table_fixed.dat" "$tmp/$nxx/MACE_RELAXATION.dat")
            if [ -n "$check" ]; then
                fail_line "$label" "[-writerelaxed no: $check]"
                mismatch=1
            fi
        fi
    fi

    # Without -relax there is no relaxed structure to suppress.
    out=$(cd "$tmp/$nxx" && PATH="$BIN:$PATH" sod_mace.sh -device cpu -cueq off \
          -writerelaxed no -force -q 2>&1)
    if [ $? -eq 0 ] || ! printf '%s' "$out" | grep -q "writerelaxed only applies with relax"; then
        fail_line "$label" "[-writerelaxed no accepted without -relax]"
        mismatch=1
    fi

    rm -rf "$tmp"
    if [ $mismatch -eq 0 ]; then
        pass_line "$label"; pass=$((pass+1))
    else
        fail=$((fail+1))
    fi
}

test_sod_mace_relaxcell() {
    local label="$1" xdir="$2" nxx="$3"
    local py="${SOD_PYTHON:-python3}"

    if [ ! -f "$xdir/INSOD" ] || [ ! -f "$xdir/SGO" ]; then
        skip_line "$label" "(missing INSOD or SGO)"
        skip=$((skip+1)); return
    fi
    if ! command -v "$py" >/dev/null 2>&1 || \
       ! "$py" -c 'import torch, ase, mace, nvalchemi' >/dev/null 2>&1; then
        skip_line "$label" "(python MACE stack not installed)"
        skip=$((skip+1)); return
    fi

    # A position-only skin cache is unsafe when the periodic cell moves: new
    # image neighbours can enter the cutoff while Cartesian atoms stay fixed.
    # The production hook factory must therefore disable the skin for
    # variable-cell work. This tiny geometry reproduces the original bug.
    local neighbor_check
    neighbor_check=$("$py" - "$ROOT" 2>/dev/null <<'PYEOF'
import pathlib
import sys

import torch
from ase import Atoms
from nvalchemi.data import AtomicData, Batch
from nvalchemi.models.base import NeighborConfig

sys.path.insert(0, str(pathlib.Path(sys.argv[1]) / "pysod"))
import mace_backend as mb

atoms = Atoms("HH", positions=[[0, 0, 0], [4, 0, 0]], cell=[10, 10, 10], pbc=True)
data = AtomicData.from_atoms(atoms, device="cpu", dtype=torch.float32)
batch = Batch.from_data_list([data], device="cpu")
hook = mb._neighbor_hook(NeighborConfig(cutoff=3.0), 0.3, variable_cell=True)
if hook.skin != 0.0:
    raise SystemExit(f"variable-cell neighbour skin is {hook.skin}, expected 0")
mb._rebuild_edges(batch, hook)
if batch.neighbor_list.shape[0] != 0:
    raise SystemExit("test geometry unexpectedly has initial neighbours")
batch.cell[0, 0, 0] = 5.0
mb._rebuild_edges(batch, hook)
if batch.neighbor_list.shape[0] != 2:
    raise SystemExit(
        f"cell-only contraction produced {batch.neighbor_list.shape[0]} directed neighbours, expected 2"
    )
PYEOF
    )
    if [ $? -ne 0 ]; then
        fail_line "$label" "[cell-aware neighbour-list regression: $neighbor_check]"
        fail=$((fail+1)); return
    fi

    local tmp; tmp=$(mktemp -d)
    cp "$xdir/INSOD" "$xdir/SGO" "$tmp/"
    mkdir -p "$tmp/$nxx"
    local n=0 c
    for c in "$xdir/$nxx"/c*/; do
        [ -d "$c" ] || continue
        cp -r "$c" "$tmp/$nxx/" || true
        # Drop any relaxed structure a previous sod_mace.sh run left in the
        # example tree: copying it in would trip the result-file protection.
        rm -f "$tmp/$nxx/$(basename "${c%/}")"/relaxed.* \
              "$tmp/$nxx/$(basename "${c%/}")"/relaxed_*
        n=$((n+1))
        [ $n -ge 2 ] && break
    done
    # This focused test intentionally copies only two configurations, so give
    # it a matching two-row ENSEMBLE. CELL is valid only when these agree.
    awk -v want="$n" '
        NR == 1 { sub(/[0-9]+ configurations/, want " configurations") }
        data && $1 ~ /^[0-9]+$/ { rows++; if (rows > want) next }
        /^# Configuration/ { data=1 }
        { print }
    ' "$xdir/$nxx/ENSEMBLE" > "$tmp/$nxx/ENSEMBLE"

    # Few steps: this checks the machinery and the direction of motion, not
    # convergence, which would be slow on CPU.
    local out rc
    out=$(cd "$tmp/$nxx" && PATH="$BIN:$PATH" sod_mace.sh -device cpu -cueq off \
          -batch 2 -relax -relaxcell -lattice cub -maxsteps 15 -force -q 2>&1)
    rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[-relaxcell error]"
        echo "$out" | tail -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local mismatch=0 check
    # CELL: one row per configuration, cubic columns "a V", no index column.
    check=$(awk -v want="$n" '
        /^[[:space:]]*(#|$)/ { next }
        { rows++; if (NF != 2) { printf "row %d has %d columns, expected 2\n", rows, NF; exit }
          if ($1 <= 0 || $2 <= 0) { printf "non-positive value in row %d\n", rows; exit } }
        END { if (rows != want) printf "expected %d rows, got %d\n", want, rows }
    ' "$tmp/$nxx/CELL" 2>/dev/null)
    if [ ! -f "$tmp/$nxx/CELL" ]; then
        fail_line "$label" "[no CELL written]"; mismatch=1
    elif [ -n "$check" ]; then
        fail_line "$label" "[CELL: $check]"; mismatch=1
    fi
    rm -f "$tmp/$nxx/ENERGIES"
    out=$(cd "$tmp/$nxx" && PATH="$BIN:$PATH" sod_mace.sh -device cpu -cueq off \
          -batch 2 -relax -relaxcell -lattice cub -maxsteps 15 -q 2>&1)
    if [ $? -eq 0 ] || ! printf '%s' "$out" | grep -q "CELL"; then
        fail_line "$label" "[rerun did not refuse existing CELL]"
        mismatch=1
    fi

    # The table must carry the variable-cell columns, and the cell must move the
    # right way: this system is under positive pressure, so it expands at P=0.
    check=$("$py" - "$tmp/$nxx" <<'PYEOF'
import sys, pathlib
from ase.io import read
level = pathlib.Path(sys.argv[1])
head, *rows = [l.split() for l in (level/"MACE_RELAXATION.dat").read_text().splitlines() if l.strip()]
for col in ("initial_volume_A3", "final_volume_A3", "final_pressure_GPa"):
    if col not in head:
        print(f"{col} missing from MACE_RELAXATION.dat header"); sys.exit()
vi, vf = head.index("initial_volume_A3") - 1, head.index("final_volume_A3") - 1
for row in rows:
    if float(row[vf]) <= float(row[vi]):
        print(f"config {row[0]}: volume did not expand at P=0 "
              f"({row[vi]} -> {row[vf]})"); sys.exit()
for d in sorted(level.glob("c*")):
    if not d.is_dir():
        continue
    a = read(d/"configuration.cif").cell.lengths()[0]
    b = read(d/"relaxed.cif").cell.lengths()[0]
    if abs(a - b) < 1e-6:
        print(f"{d.name}: relaxed.cif kept the input cell ({a})"); sys.exit()
PYEOF
)
    if [ -n "$check" ]; then
        fail_line "$label" "[$check]"; mismatch=1
    fi
    if [ -f "$tmp/$nxx/ENTHALPIES" ]; then
        fail_line "$label" "[zero-pressure MACE run wrote ENTHALPIES]"
        mismatch=1
    fi

    # Sign guard: a large target pressure must compress instead. FIRE2VariableCell
    # reads the cell force from batch.stress, so a target pressure that never
    # reaches the stored stress would silently expand here too.
    out=$(cd "$tmp/$nxx" && PATH="$BIN:$PATH" sod_mace.sh -device cpu -cueq off \
          -batch 2 -relax -relaxcell -pressure 50 -maxsteps 15 -force -q 2>&1)
    if [ $? -ne 0 ]; then
        fail_line "$label" "[-pressure error]"
        echo "$out" | tail -3 | indent
        mismatch=1
    else
        check=$(awk '
            /^#/ { for (i=1;i<=NF;i++) { if ($i=="initial_volume_A3") vi=i-1
                                         if ($i=="final_volume_A3") vf=i-1 } ; next }
            NF { if ($vf >= $vi) { printf "config %s: volume did not compress at 50 GPa (%s -> %s)\n", $1, $vi, $vf; exit } }
        ' "$tmp/$nxx/MACE_RELAXATION.dat")
        if [ -n "$check" ]; then
            fail_line "$label" "[$check]"; mismatch=1
        fi
        check=$("$py" - "$tmp/$nxx" 2>&1 <<'PYEOF'
import pathlib
import sys

GPA_PER_EV_A3 = 160.21766208
level = pathlib.Path(sys.argv[1])
energy_lines = (level / "ENERGIES").read_text().splitlines()
enthalpy_lines = (level / "ENTHALPIES").read_text().splitlines()
if "internal energy E" not in energy_lines[0]:
    print("ENERGIES header does not identify values as internal energies")
    raise SystemExit(1)
if "enthalpy H=E+PV" not in enthalpy_lines[0]:
    print("ENTHALPIES header does not identify values as enthalpies")
    raise SystemExit(1)
energies = {
    int(words[0]): float(words[1])
    for line in energy_lines
    if line.strip() and not line.lstrip().startswith("#")
    for words in [line.split()]
}
enthalpies = {
    int(words[0]): float(words[1])
    for line in enthalpy_lines
    if line.strip() and not line.lstrip().startswith("#")
    for words in [line.split()]
}
lines = [line.split() for line in (level / "MACE_RELAXATION.dat").read_text().splitlines()]
header, rows = lines[0], lines[1:]
energy_col = header.index("final_energy_eV") - 1
volume_col = header.index("final_volume_A3") - 1
for row in rows:
    index = int(row[0])
    energy = float(row[energy_col])
    expected_h = energy + 50.0 * float(row[volume_col]) / GPA_PER_EV_A3
    if abs(energies[index] - energy) > 1.0e-4:
        print(
            f"config {index}: ENERGIES {energies[index]:.8f}, "
            f"expected internal energy {energy:.8f}"
        )
        raise SystemExit(1)
    if abs(enthalpies[index] - expected_h) > 1.0e-4:
        print(
            f"config {index}: ENTHALPIES {enthalpies[index]:.8f}, "
            f"expected E+PV {expected_h:.8f}"
        )
        raise SystemExit(1)
PYEOF
        )
        rc=$?
        if [ $rc -ne 0 ] || [ -n "$check" ]; then
            fail_line "$label" "[$check]"; mismatch=1
        fi
    fi

    # MACE follows the same general policy: even valid cells for a sparse
    # subset must not become an unindexed, apparently dense CELL file.
    check=$("$py" - "$ROOT" "$tmp/$nxx" 2>&1 <<'PYEOF'
import pathlib
import sys
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(pathlib.Path(sys.argv[1]) / "pysod"))
import sod_mace

level = pathlib.Path(sys.argv[2])
result = SimpleNamespace(cell=np.eye(3))
sod_mace.write_cell_file(level, {1: result}, "cub", lambda message: None)
if (level / "CELL").exists():
    raise SystemExit("sparse MACE results wrote CELL")
PYEOF
    )
    if [ $? -ne 0 ] || ! printf '%s' "$check" | grep -q "1 of 2 ENSEMBLE configurations"; then
        fail_line "$label" "[sparse MACE results were allowed to write CELL]"
        mismatch=1
    fi

    rm -rf "$tmp"
    if [ $mismatch -eq 0 ]; then
        pass_line "$label"; pass=$((pass+1))
    else
        fail=$((fail+1))
    fi
}

test_sod_mace_settings() {
    local label="$1" xdir="$2" nxx="$3"
    local py="${SOD_PYTHON:-python3}"

    if [ ! -f "$xdir/INSOD" ] || [ ! -f "$xdir/SGO" ]; then
        skip_line "$label" "(missing INSOD or SGO)"
        skip=$((skip+1)); return
    fi
    if ! command -v "$py" >/dev/null 2>&1 || \
       ! "$py" -c 'import torch, ase, mace, nvalchemi, yaml' >/dev/null 2>&1; then
        skip_line "$label" "(python MACE stack not installed)"
        skip=$((skip+1)); return
    fi

    local tmp; tmp=$(mktemp -d)
    cp "$xdir/INSOD" "$xdir/SGO" "$tmp/"
    mkdir -p "$tmp/$nxx"
    local n=0 c
    for c in "$xdir/$nxx"/c*/; do
        [ -d "$c" ] || continue
        cp -r "$c" "$tmp/$nxx/" || true
        # Drop any relaxed structure a previous sod_mace.sh run left in the
        # example tree: copying it in would trip the result-file protection.
        rm -f "$tmp/$nxx/$(basename "${c%/}")"/relaxed.* \
              "$tmp/$nxx/$(basename "${c%/}")"/relaxed_*
        n=$((n+1))
        [ $n -ge 2 ] && break
    done

    # Values chosen so the run is only explicable by the file having been read:
    # they are echoed verbatim into the ENERGIES provenance header.
    cat > "$tmp/$nxx/mace_settings.yaml" <<'YAML'
device: cpu
cueq: off
batch: 2
relax: true
batchmode: refill
fmax: 0.06
maxsteps: 20
q: true
YAML

    local mismatch=0 out rc
    out=$(cd "$tmp/$nxx" && PATH="$BIN:$PATH" sod_mace.sh 2>&1)
    rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[sod_mace.sh error with settings file]"
        echo "$out" | tail -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local header; header=$(head -1 "$tmp/$nxx/ENERGIES" 2>/dev/null)
    case "$header" in
        *"refill batching"*"fmax 0.06"*"max 20 steps"*) ;;
        *)  fail_line "$label" "[settings not applied]"
            printf '%s\n' "$header" | indent
            mismatch=1 ;;
    esac

    # An unknown key must be rejected, not silently ignored.
    echo "batchsize: 4" >> "$tmp/$nxx/mace_settings.yaml"
    out=$(cd "$tmp/$nxx" && PATH="$BIN:$PATH" sod_mace.sh -force 2>&1)
    if [ $? -eq 0 ] || ! printf '%s' "$out" | grep -q "unknown option"; then
        fail_line "$label" "[unknown settings key not rejected]"
        mismatch=1
    fi

    rm -rf "$tmp"
    if [ $mismatch -eq 0 ]; then
        pass_line "$label"; pass=$((pass+1))
    else
        fail=$((fail+1))
    fi
}

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

# ── test_cpmesod ───────────────────────────────────────────────────────────────
# $1 = display label
# $2 = main example directory (contains INSOD, SGO, n00-n03, n21-n24)
# $3 = target level
# Reference file must already be committed in $2/cpme_test_ref/ENERGIES
test_cpmesod() {
    local label="$1" maindir="$2" tlvl="$3"

    if [ ! -f "$maindir/INSOD" ] || [ ! -f "$maindir/SGO" ]; then
        skip_line "$label" "(missing INSOD or SGO)"
        skip=$((skip+1)); return
    fi

    # Low-side levels: 0-3; high-side levels: 21-24

    # Check that all required level directories and data files exist
    local missing=0
    for i in 0 1 2 3 21 22 23 24; do
        local ndir; ndir=$(sod_level_dir_by_number "$maindir" $i)
        if [ -z "$ndir" ] || [ ! -f "$maindir/$ndir/ENSEMBLE" ] || [ ! -f "$maindir/$ndir/ENERGIES" ]; then
            fail_line "$label" "[level $i: ENSEMBLE or ENERGIES missing]"
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
        local ndir; ndir=$(sod_level_dir_by_number "$maindir" $i)
        mkdir -p "$tmp/$ndir"
        cp "$maindir/$ndir/ENSEMBLE" "$maindir/$ndir/ENERGIES" "$tmp/$ndir/"
    done

    # Run combsod to generate EQMATRIX
    local out; out=$(cd "$tmp" && PATH="$BIN:$PATH" combsod 2>&1)
    local rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[combsod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Run cpmesod
    out=$(cd "$tmp" && PATH="$BIN:$PATH" cpmesod 2>&1)
    rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[cpmesod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Check CPMEh ENERGIES from target level
    local tdir; tdir=$(sod_level_dir_by_number "$maindir" $tlvl)
    local ref="$maindir/cpme_test_ref/ENERGIES"
    local gen="$tmp/$tdir/CPMEh/ENERGIES"
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

# ── test_cpmesod_example17 ─────────────────────────────────────────────────────
# example17: n00–n03 (low, order_cap=3) + n24–n27 (high), target n04/CPMEh.
# cpme.model is copied so cpmesod sees cpme_ch=2/order_cap=3/n_calib=0/alpha=2/eta=0.
# n04 is the prediction target, not a training level.
# n04/ENSEMBLE and n04/ENERGIES are deliberately included in the tmp dir to
# verify that cpmesod does NOT use them as training data (the target-level guard).
# Without the guard, cpmesod would crash with a "training data conflict" error.
# n_calib=0 → manual eps mode (eps=1, alpha=2, eta=0): "epsilon read/fitted from cpme.model."
# $1 = display label
# $2 = main example directory
test_cpmesod_example17() {
    local label="$1" maindir="$2"

    if [ ! -f "$maindir/INSOD" ] || [ ! -f "$maindir/SGO" ] || [ ! -f "$maindir/cpme.model" ]; then
        skip_line "$label" "(missing INSOD, SGO or cpme.model)"
        skip=$((skip+1)); return
    fi

    local missing=0
    for i in 0 1 2 3 4 24 25 26 27; do
        local ndir; ndir=$(sod_level_dir_by_number "$maindir" $i)
        if [ -z "$ndir" ] || [ ! -f "$maindir/$ndir/ENSEMBLE" ] || [ ! -f "$maindir/$ndir/ENERGIES" ]; then
            fail_line "$label" "[level $i: ENSEMBLE or ENERGIES missing]"
            missing=1; break
        fi
    done
    [ $missing -ne 0 ] && { fail=$((fail+1)); return; }

    local tmp; tmp=$(mktemp -d)
    cp "$maindir/INSOD" "$maindir/SGO" "$maindir/cpme.model" "$tmp/"

    # Include n04 so the target-level guard is exercised: cpmesod must not use
    # n04/ENERGIES as a training reference even though the file is present.
    for i in 0 1 2 3 4 24 25 26 27; do
        local ndir; ndir=$(sod_level_dir_by_number "$maindir" $i)
        mkdir -p "$tmp/$ndir"
        cp "$maindir/$ndir/ENSEMBLE" "$maindir/$ndir/ENERGIES" "$tmp/$ndir/"
    done

    local out; out=$(cd "$tmp" && PATH="$BIN:$PATH" combsod 2>&1)
    local rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[combsod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    out=$(cd "$tmp" && PATH="$BIN:$PATH" cpmesod 2>&1)
    rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[cpmesod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local tdir; tdir=$(sod_level_dir_by_number "$maindir" 4)
    local ref="$maindir/cpme_test_ref/ENERGIES"
    local gen="$tmp/$tdir/CPMEh/ENERGIES"
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

# ── test_cpme_delta ────────────────────────────────────────────────────────────
# $1 = display label
# $2 = main example directory (contains INSOD, SGO, n*/ENSEMBLE+ENERGIES)
# $3 = space-separated level list to copy (low + high side)
# Builds the same CPME model mcsod would use (combsod + cpmesod), then runs the
# test_cpme_delta driver, which checks the incremental cpme_evaluate_swap_delta
# against the full cpme_evaluate_configuration over a level sweep and a long
# swap chain. The driver exits non-zero on any mismatch.
test_cpme_delta_case() {
    local label="$1" maindir="$2" levels="$3"

    if [ ! -f "$maindir/INSOD" ] || [ ! -f "$maindir/SGO" ]; then
        skip_line "$label" "(missing INSOD or SGO)"
        skip=$((skip+1)); return
    fi
    if [ ! -x "$BIN/test_cpme_delta" ]; then
        skip_line "$label" "(test_cpme_delta not built — run 'make testbin')"
        skip=$((skip+1)); return
    fi

    local i ndir missing=0
    for i in $levels; do
        ndir=$(sod_level_dir_by_number "$maindir" "$i")
        if [ -z "$ndir" ] || [ ! -f "$maindir/$ndir/ENSEMBLE" ] || [ ! -f "$maindir/$ndir/ENERGIES" ]; then
            fail_line "$label" "[level $i: ENSEMBLE or ENERGIES missing]"
            missing=1; break
        fi
    done
    [ $missing -ne 0 ] && { fail=$((fail+1)); return; }

    local tmp; tmp=$(mktemp -d)
    cp "$maindir/INSOD" "$maindir/SGO" "$tmp/"
    [ -f "$maindir/cpme.model" ] && cp "$maindir/cpme.model" "$tmp/"
    for i in $levels; do
        ndir=$(sod_level_dir_by_number "$maindir" "$i")
        mkdir -p "$tmp/$ndir"
        cp "$maindir/$ndir/ENSEMBLE" "$maindir/$ndir/ENERGIES" "$tmp/$ndir/"
    done

    local out rc
    out=$(cd "$tmp" && PATH="$BIN:$PATH" combsod 2>&1); rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[combsod error]"; echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    out=$(cd "$tmp" && PATH="$BIN:$PATH" cpmesod 2>&1); rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[cpmesod error]"; echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    out=$(cd "$tmp" && PATH="$BIN:$PATH" test_cpme_delta 2>&1); rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[delta != full]"
        echo "$out" | grep -iE "mismatch|drift|max |error" | head -8 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    rm -rf "$tmp"
    pass_line "$label"; pass=$((pass+1))
}

# ── test_mcsod ────────────────────────────────────────────────────────────────
# $1 = display label
# $2 = main example directory (contains INSOD, SGO, INMC, n00-n03, n21-n24)
# $3 = target level (where OUTMC is written)
# Reference file must already be committed in $2/cpme_ref/OUTMC
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
        local ndir; ndir=$(sod_level_dir_by_number "$maindir" $i)
        if [ -z "$ndir" ] || [ ! -f "$maindir/$ndir/ENSEMBLE" ] || [ ! -f "$maindir/$ndir/ENERGIES" ]; then
            fail_line "$label" "[level $i: ENSEMBLE or ENERGIES missing]"
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
        local ndir; ndir=$(sod_level_dir_by_number "$maindir" $i)
        mkdir -p "$tmp/$ndir"
        cp "$maindir/$ndir/ENSEMBLE" "$maindir/$ndir/ENERGIES" "$tmp/$ndir/"
    done

    # Run combsod then cpmesod to generate Hamiltonian (needed for mcsod)
    local out; out=$(cd "$tmp" && PATH="$BIN:$PATH" combsod 2>&1)
    local rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[combsod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    out=$(cd "$tmp" && PATH="$BIN:$PATH" cpmesod 2>&1)
    rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[cpmesod error]"
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
    local tdir; tdir=$(sod_level_dir_by_number "$maindir" $tlvl)
    local ref="$maindir/cpme_test_ref/OUTMC"
    local gen="$tmp/$tdir/MCT_300K/CPMEh/OUTMC"
    if [ ! -f "$ref" ]; then
        fail_line "$label" "[OUTMC has no reference]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    elif ! mc_outmc_match "$ref" "$gen" 0.02; then
        fail_line "$label" "[OUTMC differs beyond tolerance]"
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
# Reference file must already be committed in $2/cpme_test_ref/thermodynamics.dat
test_mcstat() {
    local label="$1" maindir="$2" tlvl="$3"

    if [ ! -f "$maindir/INSOD" ] || [ ! -f "$maindir/SGO" ] || [ ! -f "$maindir/INMC" ]; then
        skip_line "$label" "(missing INSOD, SGO, or INMC)"
        skip=$((skip+1)); return
    fi

    # Low-side levels: 0-3; high-side levels: 21-24
    local missing=0
    for i in 0 1 2 3 21 22 23 24; do
        local ndir; ndir=$(sod_level_dir_by_number "$maindir" $i)
        if [ -z "$ndir" ] || [ ! -f "$maindir/$ndir/ENSEMBLE" ] || [ ! -f "$maindir/$ndir/ENERGIES" ]; then
            fail_line "$label" "[level $i: ENSEMBLE or ENERGIES missing]"
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
        local ndir; ndir=$(sod_level_dir_by_number "$maindir" $i)
        mkdir -p "$tmp/$ndir"
        cp "$maindir/$ndir/ENSEMBLE" "$maindir/$ndir/ENERGIES" "$tmp/$ndir/"
    done

    # Run combsod then cpmesod to generate the Hamiltonian (needed for mcsod)
    local out; out=$(cd "$tmp" && PATH="$BIN:$PATH" combsod 2>&1)
    local rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[combsod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    out=$(cd "$tmp" && PATH="$BIN:$PATH" cpmesod 2>&1)
    rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[cpmesod error]"
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

    # Thermodynamic integration over the sampled temperatures.
    # Layout is sampling-first: MC output lives in nXX/MCT_*K/CPMEx/, so mcstatsod
    # runs from nXX/ and writes thermodynamics.dat there.
    local tdir; tdir=$(sod_level_dir_by_number "$maindir" $tlvl)
    local nxxdir="$tmp/$tdir"
    if ! ls -d "$nxxdir"/MCT_*K/ >/dev/null 2>&1; then
        fail_line "$label" "[$tdir/MCT_*K not generated]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    out=$(cd "$nxxdir" && PATH="$BIN:$PATH" sod_mcstat.sh 2>&1)
    rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[mcstatsod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Check thermodynamics.dat against the committed reference
    local ref="$maindir/cpme_test_ref/thermodynamics.dat"
    local gen="$nxxdir/thermodynamics.dat"
    if [ ! -f "$ref" ]; then
        fail_line "$label" "[thermodynamics.dat has no reference]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    elif ! thermo_match "$ref" "$gen"; then
        fail_line "$label" "[thermodynamics.dat differs beyond tolerance]"
        diff "$ref" "$gen" | head -6 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    rm -rf "$tmp"
    pass_line "$label"; pass=$((pass+1))
}

# ── test_mcstat_vs_enum ───────────────────────────────────────────────────────
# Cross-check: thermodynamic integration (mcstatsod) vs the exact full
# enumeration (statsod) on the *same* CPMEh Hamiltonian. Both paths use the CPMEh
# energies that cpmesod predicts for every configuration of the target level, so
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
        local ndir; ndir=$(sod_level_dir_by_number "$maindir" $i)
        if [ -z "$ndir" ] || [ ! -f "$maindir/$ndir/ENSEMBLE" ] || [ ! -f "$maindir/$ndir/ENERGIES" ]; then
            fail_line "$label" "[level $i: ENSEMBLE or ENERGIES missing]"
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
        local ndir; ndir=$(sod_level_dir_by_number "$maindir" $i)
        mkdir -p "$tmp/$ndir"
        cp "$maindir/$ndir/ENSEMBLE" "$maindir/$ndir/ENERGIES" "$tmp/$ndir/"
    done

    # combsod (full enumeration of target level) + cpmesod (CPMEh energies)
    local out rc
    out=$(cd "$tmp" && PATH="$BIN:$PATH" combsod 2>&1); rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[combsod error]"; echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    out=$(cd "$tmp" && PATH="$BIN:$PATH" cpmesod 2>&1); rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[cpmesod error]"; echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local tdir; tdir=$(sod_level_dir_by_number "$maindir" $tlvl)
    local nxxdir="$tmp/$tdir"
    local cpmedir="$nxxdir/CPMEh"   # cpmesod enumeration energies (nXX/CPMEh/ENERGIES)
    if [ ! -f "$cpmedir/ENERGIES" ] || [ ! -f "$nxxdir/ENSEMBLE" ]; then
        fail_line "$label" "[$tdir enumeration/ENERGIES not generated]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Exact path: statsod over the full enumeration with the CPMEh energies.
    # cpmesod writes the indexed two-column "m  E" form statsod reads, so the
    # file goes straight across (before 0.91 it was a bare energy column and
    # this test had to prepend the index itself).
    mkdir -p "$tmp/enum"
    cp "$nxxdir/ENSEMBLE" "$tmp/enum/ENSEMBLE"
    cp "$tmp/TEMPERATURES"      "$tmp/enum/TEMPERATURES"
    cp "$cpmedir/ENERGIES"      "$tmp/enum/ENERGIES"
    out=$(cd "$tmp/enum" && PATH="$BIN:$PATH" statsod 2>&1); rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[statsod error]"; echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Approximate path: mcsod MC sampling (same CPMEh) + mcstatsod TI.
    # MC output is sampling-first (nXX/MCT_*K/CPMEh/), so mcstatsod runs from nXX/.
    out=$(cd "$tmp" && PATH="$BIN:$PATH" sod_mc.sh 2>&1); rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[mcsod error]"; echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    out=$(cd "$nxxdir" && PATH="$BIN:$PATH" sod_mcstat.sh 2>&1); rc=$?
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
    ' "$tmp/enum/thermodynamics.dat" "$nxxdir/thermodynamics.dat")
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
    cp "$maindir/INSOD" "$maindir/SGO" "$tmp/"

    # Copy n*/ENSEMBLE+ENERGIES+DATA+SPECTRA for the required range
    local missing=0
    for ((i=nsubsmin; i<=nsubsmax; i++)); do
        local ndir; ndir=$(sod_level_dir_by_number "$maindir" $i)
        if [ -z "$ndir" ] || [ ! -f "$maindir/$ndir/ENSEMBLE" ] || [ ! -f "$maindir/$ndir/ENERGIES" ]; then
            fail_line "$label" "[level $i: ENSEMBLE or ENERGIES missing]"
            missing=1; break
        fi
        mkdir -p "$tmp/$ndir"
        cp "$maindir/$ndir/ENSEMBLE" "$maindir/$ndir/ENERGIES" "$tmp/$ndir/"
        [ -f "$maindir/$ndir/DATA"    ] && cp "$maindir/$ndir/DATA"    "$tmp/$ndir/"
        [ -f "$maindir/$ndir/SPECTRA" ] && cp "$maindir/$ndir/SPECTRA" "$tmp/$ndir/"
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

# ── test_sqssod_general ───────────────────────────────────────────────────────
# Smoke-test generalized SQS paths without committing bulky reference outputs.
# $1 = display label
# $2 = main example directory (contains INSOD, EQMATRIX, supercell.cif)
# $3 = nXX directory name
# $4 = extended grep pattern expected in OUTSQS
# $5 = extended grep pattern expected in SQS_CORRELATIONS
test_sqssod_general() {
    local label="$1" maindir="$2" nxx="$3" out_pat="$4" corr_pat="$5"
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

    local tmp; tmp=$(mktemp -d)
    cp "$maindir/INSOD" "$maindir/EQMATRIX" "$maindir/supercell.cif" "$tmp/"
    mkdir -p "$tmp/$nxx"
    cp "$ndir/ENSEMBLE" "$tmp/$nxx/"
    cat > "$tmp/INSQS" <<'EOF'
# Maximum cluster order
2

# Cutoff radii (Angstroms) for orders 2..MaxOrder
8.0

# Weights for orders 2..MaxOrder
1.0

# omega and eps_tol for van de Walle scoring
10  1.0E-6

# n_top_sqs: number of top configurations listed in OUTSQS (0 = rank and list all)
10
EOF

    local out; out=$(cd "$tmp" && PATH="$BIN:$PATH" sqssod "$nxx" 2>&1)
    local rc=$?
    if [ $rc -ne 0 ]; then
        fail_line "$label" "[sqssod error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    if [ ! -s "$tmp/$nxx/OUTSQS" ] || [ ! -s "$tmp/$nxx/SQS_CORRELATIONS" ]; then
        fail_line "$label" "[missing generalized SQS outputs]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    if ! grep -Eq "$out_pat" "$tmp/$nxx/OUTSQS"; then
        fail_line "$label" "[OUTSQS missing expected family columns]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    if ! grep -Eq "$corr_pat" "$tmp/$nxx/SQS_CORRELATIONS"; then
        fail_line "$label" "[SQS_CORRELATIONS missing expected species]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    rm -rf "$tmp"
    pass_line "$label"; pass=$((pass+1))
}

# ── test_species_cap ─────────────────────────────────────────────────────────
# Validates the shared INSOD species-per-target cap: up to 5 substituent
# species (senary disorder, 6 species total) parse and flow through
# randomsod -> sqssod; a 6th substituent is rejected with a clear error
# instead of being silently truncated; and combsod's full-enumeration path
# (structurally limited to 3 substituents) rejects 4+ with a clear error.
# $1 = display label   $2 = main example dir (INSOD, SGO, EQMATRIX, supercell.cif)
test_species_cap() {
    local label="$1" maindir="$2"

    for f in SGO EQMATRIX supercell.cif; do
        if [ ! -f "$maindir/$f" ]; then
            skip_line "$label" "(missing $f)"
            skip=$((skip+1)); return
        fi
    done

    local tmp; tmp=$(mktemp -d)
    cp "$maindir/SGO" "$maindir/EQMATRIX" "$maindir/supercell.cif" "$tmp/"
    cat > "$tmp/INSOD" <<'EOF'
# Title
Senary species-cap smoke test (CoCrFeMnCu on Ni, fcc 2x2x2)

# a, b, c, alpha, beta, gamma
3.5400  3.5400  3.5400  90.000  90.000  90.000

# nsp: Number of species in the parent structure
1

# symbol(1:nsp): Atom symbols
Ni

# natsp0(1:nsp): Number of atoms per species (asymmetric unit is sufficient)
1

# coords0: Fractional coordinates
0.0  0.0  0.0

# na, nb, nc: Supercell multipliers along a, b, c
2 2 2

# sptarget: Index of species to substitute
1

# nsubs: Substitution counts for Co, Cr, Fe, Mn, Cu; remaining sites are Ni.
6 5 5 5 5

# newsymbol: new species followed by remaining species
Co Cr Fe Mn Cu Ni

# FILER: -1 none | 0 CIF | 1 GULP | 2 LAMMPS | 11 VASP | 12 CASTEP | 13 QE
0
EOF

    # 5 substituents (6 species total): randomsod -sym on must succeed and
    # produce the expected per-species ENSEMBLE header.
    local out
    out=$(cd "$tmp" && PATH="$BIN:$PATH" randomsod -nconf 200 -sym on -seed 20260710 2>&1)
    if [ $? -ne 0 ]; then
        fail_line "$label" "[randomsod rejected a valid 5-substituent target]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    local ens="$tmp/n06_05_05_05_05/random/ENSEMBLE"
    if ! grep -q "Co.*Cr.*Fe.*Mn.*Cu.*Ni" "$ens" 2>/dev/null; then
        fail_line "$label" "[ENSEMBLE missing expected 6-species target line]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # sqssod on that ensemble: 4 fcc pair orbits within the default cutoff,
    # each with 6x6=36 species channels.
    cat > "$tmp/INSQS" <<'EOF'
# Maximum cluster order
2

# Cutoff radii (Angstroms) for orders 2..MaxOrder
6.0

# Weights for orders 2..MaxOrder
1.0

# omega and eps_tol for van de Walle scoring
10  1.0E-6

# n_top_sqs: number of top configurations listed in OUTSQS (0 = rank and list all)
10
EOF
    out=$(cd "$tmp" && PATH="$BIN:$PATH" sqssod n06_05_05_05_05/random 2>&1)
    if [ $? -ne 0 ]; then
        fail_line "$label" "[sqssod error on 6-species ensemble]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    if ! grep -q "n_channels=144" "$tmp/n06_05_05_05_05/random/OUTSQS"; then
        fail_line "$label" "[expected n_channels=144 (4 orbits x 6x6 species)]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # 6 substituents (7 species) on the same target: must be rejected with a
    # clear error, not silently truncated to 5 (dropping the 6th species).
    sed 's/^6 5 5 5 5$/6 5 5 5 5 5/; s/^Co Cr Fe Mn Cu Ni$/Co Cr Fe Mn Cu Zn Ni/' \
        "$tmp/INSOD" > "$tmp/INSOD.overflow"
    mv "$tmp/INSOD.overflow" "$tmp/INSOD"
    out=$(cd "$tmp" && PATH="$BIN:$PATH" randomsod -nconf 10 -sym off -seed 1 2>&1)
    if [ $? -eq 0 ] || ! echo "$out" | grep -q "at most 5 substituent species"; then
        fail_line "$label" "[6-substituent target not rejected with expected error]"
        echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # combsod: 5 substituents (valid for the shared parser) must still be
    # rejected because full enumeration only supports up to 3.
    sed 's/^6 5 5 5 5 5$/6 5 5 5 5/; s/^Co Cr Fe Mn Cu Zn Ni$/Co Cr Fe Mn Cu Ni/' \
        "$tmp/INSOD" > "$tmp/INSOD.quinary"
    mv "$tmp/INSOD.quinary" "$tmp/INSOD"
    out=$(cd "$tmp" && PATH="$BIN:$PATH" combsod 2>&1)
    if [ $? -eq 0 ] || ! echo "$out" | grep -q "at most 3 substituent species"; then
        fail_line "$label" "[combsod did not reject a 5-substituent target]"
        echo "$out" | head -3 | indent
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

# ── test_randomsod ────────────────────────────────────────────────────────────
# Smoke test for the uniform random sampler: run randomsod with a fixed seed,
# check it writes nXX/random/ENSEMBLE, that sum_degeneracies == nconfigs (visit
# counts over the draws), and that a fixed seed is reproducible.
# $1 = display label  $2 = main example dir (needs INSOD, SGO, EQMATRIX)
test_randomsod() {
    local label="$1" maindir="$2"

    if [ ! -f "$maindir/INSOD" ] || [ ! -f "$maindir/SGO" ] || [ ! -f "$maindir/EQMATRIX" ]; then
        skip_line "$label" "(missing INSOD, SGO, or EQMATRIX)"
        skip=$((skip+1)); return
    fi

    local nc=5000
    local tmp; tmp=$(mktemp -d)
    cp "$maindir/INSOD" "$maindir/SGO" "$maindir/EQMATRIX" "$tmp/"
    local out; out=$(cd "$tmp" && PATH="$BIN:$PATH" randomsod -nconf $nc -sym on -seed 12345 2>&1)
    if [ $? -ne 0 ]; then
        fail_line "$label" "[randomsod error]"; echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local ens; ens=$(ls "$tmp"/n*/random/ENSEMBLE 2>/dev/null | head -1)
    if [ -z "$ens" ]; then
        fail_line "$label" "[no nXX/random/ENSEMBLE produced]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    local sumdeg; sumdeg=$(awk -F'sum_degeneracies = ' 'NR==1{print $2+0; exit}' "$ens")
    if [ "$sumdeg" != "$nc" ]; then
        fail_line "$label" "[sum_degeneracies=$sumdeg != nconfigs=$nc]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Reproducibility: the same fixed seed must reproduce the ENSEMBLE.
    local tmp2; tmp2=$(mktemp -d)
    cp "$maindir/INSOD" "$maindir/SGO" "$maindir/EQMATRIX" "$tmp2/"
    (cd "$tmp2" && PATH="$BIN:$PATH" randomsod -nconf $nc -sym on -seed 12345 >/dev/null 2>&1)
    local ens2; ens2=$(ls "$tmp2"/n*/random/ENSEMBLE 2>/dev/null | head -1)
    if ! diff -q "$ens" "$ens2" >/dev/null 2>&1; then
        fail_line "$label" "[not reproducible for fixed seed]"
        fail=$((fail+1)); rm -rf "$tmp" "$tmp2"; return
    fi

    rm -rf "$tmp" "$tmp2"
    pass_line "$label"; pass=$((pass+1))
}

# ── test_randomsod_general ────────────────────────────────────────────────────
# Multinary/multisite validation. combsod enumerates the exact orbit set; with a
# fixed seed and enough draws randomsod -sym on must fold to the same NUMBER of
# distinct orbits (a representative-independent invariant), with sum_degeneracies
# == nconfigs. Also checks -sym off writes one row per draw with the expected
# per-species columns.
# $1 = display label   $2 = main example dir (INSOD, SGO, EQMATRIX)
test_randomsod_general() {
    local label="$1" maindir="$2"

    if [ ! -f "$maindir/INSOD" ] || [ ! -f "$maindir/SGO" ] || [ ! -f "$maindir/EQMATRIX" ]; then
        skip_line "$label" "(missing INSOD, SGO, or EQMATRIX)"
        skip=$((skip+1)); return
    fi

    local tmp; tmp=$(mktemp -d)
    cp "$maindir/INSOD" "$maindir/SGO" "$tmp/"

    # Reference orbit count from combsod's exact enumeration.
    (cd "$tmp" && PATH="$BIN:$PATH" combsod >/dev/null 2>&1)
    local cens; cens=$(ls "$tmp"/n*/ENSEMBLE 2>/dev/null | grep -v /random/ | head -1)
    if [ -z "$cens" ]; then
        fail_line "$label" "[combsod produced no ENSEMBLE]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    local ncomb; ncomb=$(awk -F'ensemble: ' 'NR==1{split($2,a," "); print a[1]; exit}' "$cens")

    # -sym on: enough draws to sample every orbit, then compare distinct count.
    local nc=20000
    local out; out=$(cd "$tmp" && PATH="$BIN:$PATH" randomsod -nconf $nc -sym on -seed 4242 2>&1)
    if [ $? -ne 0 ]; then
        fail_line "$label" "[randomsod -sym on error]"; echo "$out" | head -3 | indent
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    local ens; ens=$(ls "$tmp"/n*/random/ENSEMBLE 2>/dev/null | head -1)
    local ndist; ndist=$(awk -F'ensemble: ' 'NR==1{split($2,a," "); print a[1]; exit}' "$ens")
    if [ "$ndist" != "$ncomb" ]; then
        fail_line "$label" "[distinct orbits=$ndist != combsod=$ncomb]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi
    local sumdeg; sumdeg=$(awk -F'sum_degeneracies = ' 'NR==1{print $2+0; exit}' "$ens")
    if [ "$sumdeg" != "$nc" ]; then
        fail_line "$label" "[sym-on sum_degeneracies=$sumdeg != nconfigs=$nc]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # Byte-for-byte representative parity: colored canonicalization must pick the
    # same orbit representative as combsod, so the sets of position columns
    # (fields 3..NF of each data row) must be identical between the two ENSEMBLEs.
    local reps_c reps_r
    reps_c=$(awk '$1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/{$1="";$2="";print}' "$cens" | sort)
    reps_r=$(awk '$1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/{$1="";$2="";print}' "$ens"  | sort)
    if [ "$reps_c" != "$reps_r" ]; then
        fail_line "$label" "[representatives differ from combsod]"
        fail=$((fail+1)); rm -rf "$tmp"; return
    fi

    # -sym off: one row per draw, sum_degeneracies == nconfigs.
    local tmp2; tmp2=$(mktemp -d)
    cp "$maindir/INSOD" "$maindir/SGO" "$tmp2/"
    (cd "$tmp2" && PATH="$BIN:$PATH" randomsod -nconf 100 -sym off -seed 4242 >/dev/null 2>&1)
    local ensoff; ensoff=$(ls "$tmp2"/n*/random/ENSEMBLE 2>/dev/null | head -1)
    # Data rows have numeric config index AND numeric degeneracy; target header
    # lines (e.g. "8 La sites -> ...") have a non-numeric second field.
    local nrows; nrows=$(awk '$1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/{c++} END{print c+0}' "$ensoff")
    local sumoff; sumoff=$(awk -F'sum_degeneracies = ' 'NR==1{print $2+0; exit}' "$ensoff")
    if [ "$nrows" != "100" ] || [ "$sumoff" != "100" ]; then
        fail_line "$label" "[sym-off rows=$nrows sum=$sumoff, expected 100/100]"
        fail=$((fail+1)); rm -rf "$tmp" "$tmp2"; return
    fi

    rm -rf "$tmp" "$tmp2"
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
test_genersod_multitarget "example13 (3 targets, vacancy on target 3)" "$EX/example13" \
    "n01_01_01/c1/configuration.cif" "$EX/example13/genersod_ref/configuration.cif"
test_genersod_molecule_second_target "example08 (@MA on the second target, LAMMPS)" "$EX/example08"
test_genersod_multinary_mol_vac "example20 (@ and % on two multi-nary targets)" "$EX/example20"
test_eqmatrix_only "example18 (-eqmatrix-only on a non-enumerable cell)" "$EX/example18"
test_gener_best "example16 (-choose bestSQS picks the OUTSQS rank 1)" "$EX/example16" "n08"
test_sqs_no_ensemble_scored "example16 (sod_sqs.sh fails when it scores nothing)" "$EX/example16" "n08"
test_gener_choose_random "example20 (-choose random N draws without replacement)" "$EX/example20"
test_wrapper_filer_tail "sod_comb.sh ignores trailing INSOD comment" "comb" "$EX/example01/FILER1_gulp" "n04/c01/input.gin"
test_wrapper_filer_tail "sod_gener.sh ignores trailing INSOD comment" "gener" "$EX/example01/FILER1_gulp" "job_sender"
test_wrapper_model_flag_error "sod_cpme.sh rejects missing -model filename" "$BIN/sod_cpme.sh" "$EX/example15"
test_wrapper_model_flag_error "sod_mc.sh rejects missing -model filename" "$BIN/sod_mc.sh" "$EX/example15"

echo ""
printf "CELL extractors  (complete-ensemble policy)\n"
printf "%s\n" "-------------------------------------------"
test_cell_completeness "VASP/GULP CELL requires complete ENSEMBLE"
test_energy_enthalpy_extractors "VASP/GULP ENERGIES and ENTHALPIES"
test_ensemble_v3_only "ENSEMBLE version 3 required everywhere"
test_mace_result_protection "sod_mace result files protected before a run"

# ── statsod tests (canonical) ────────────────────────────────────────────────

echo ""
printf "statsod  (canonical statistics)\n"
printf "%s\n" "-------------------------------------------"

test_stat "example05/n02 (energy+DATA+SPECTRA)" "$EX/example05/n02"
test_stat "example16/n08 (8043 configs, T=0/1000/1e6 K)" "$EX/example16/n08"

# ── cpmesod tests (Constrained Periodic Motif Expansion) ────────────────────────────────────

echo ""
printf "cpmesod  (Constrained Periodic Motif Expansion)\n"
printf "%s\n" "-------------------------------------------"

test_cpmesod "example15/cpmesod (n00-n03 low, n21-n24 high, target=n12)" "$EX/example15" 12
test_cpmesod_example17 "example17/cpmesod (n00-n03 low, n24-n27 high, target=n04, CPMEh)" "$EX/example17"

# Incremental swap evaluator vs full recompute (delta == full(new) - full(old))
test_cpme_delta_case "example15/cpme_delta (swap evaluator vs full)" "$EX/example15" "0 1 2 3 21 22 23 24"
test_cpme_delta_case "example17/cpme_delta (swap evaluator vs full)" "$EX/example17" "0 1 2 3 4 24 25 26 27"

# ── randomsod tests (uniform random sampling) ─────────────────────────────────

echo ""
printf "randomsod  (uniform random sampling)\n"
printf "%s\n" "-------------------------------------------"

test_randomsod "example15/randomsod (n12, 5000 draws, symmetry on)" "$EX/example15"
test_randomsod_general "example14/randomsod (multisite+multinary vs combsod)" "$EX/example14"
test_randomsod_general "example09/randomsod (multisite vs combsod)" "$EX/example09"

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
test_sqssod "example12/sqssod multitarget (n02_02)" "$EX/example12" "n02_02"
test_sqssod_general "example10/sqssod multinary target" "$EX/example10" "n04_04" "E11" "Zr|Nb"
test_sqssod_general "example14/sqssod multinary multitarget" "$EX/example14" "n01_01_01" "E11.*E12.*E22" "Ba|Mn"
test_species_cap "example19 geometry/species cap (5 ok, 6 rejected, combsod<=3)" "$EX/example19"

# ── gqssod tests (Generalized Quasirandom Structures) ────────────────────────

echo ""
printf "gqssod  (Generalized Quasirandom Structures)\n"
printf "%s\n" "-------------------------------------------"

test_gqssod "example16/gqssod (n08, T=0/1000/1e6 K)" "$EX/example16" "n08"

# ── sod_mace tests (MACE machine-learning potential) ─────────────────────────

echo ""
printf "sod_mace  (MACE machine-learning potential)\n"
printf "%s\n" "-------------------------------------------"

test_mace_private_api "nvalchemi private APIs still present"
test_sod_mace "example01/FILER0_mace (n04, CIF input, 4 configs)" "$EX/example01/FILER0_mace" "n04"
test_sod_mace_relax "example01/FILER0_mace relax fixed vs refill parity" "$EX/example01/FILER0_mace" "n04"
test_sod_mace_relaxcell "example01/FILER0_mace variable cell + CELL + pressure sign" "$EX/example01/FILER0_mace" "n04"
test_sod_mace_settings "example01/FILER0_mace mace_settings.yaml" "$EX/example01/FILER0_mace" "n04"

# ── summary ───────────────────────────────────────────────────────────────────

echo ""
echo "Results: $pass passed, $fail failed, $skip skipped"
[ $fail -eq 0 ]
