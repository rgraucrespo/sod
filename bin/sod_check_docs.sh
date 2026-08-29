#!/usr/bin/env bash
# sod_check_docs.sh — documentation consistency checks
#
# Usage: ./bin/sod_check_docs.sh
#
# README.md is the offline manual; docs/ is the same material restructured for
# readthedocs.  Both are edited by hand, so they drift.  This script checks the
# handful of facts that actually go stale when one copy is updated and the
# other is not.  Exit status is 0 if all checks pass, 1 otherwise.
#
#   version      the release string agrees across README.md, docs/conf.py,
#                VERSIONS.md and every src/*.f90 startup banner
#   examples     every examples/exampleNN/ directory is named in README.md
#                (docs/examples.rst is a curated selection, not a catalogue,
#                so it is deliberately not checked)
#   scripts      every user-facing bin/sod_*.sh is named in README.md and
#                somewhere under docs/
#   readme-toc   every in-page link in README.md resolves to a heading
#   sphinx       docs/ builds with warnings as errors (SKIP without sphinx)
#
# In the development repository, AGENTS.md carries the section-to-page map
# between README.md and docs/.

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT" || exit 1

pass=0; fail=0; skip=0
label_fmt="%-14s"
pass_line() { printf "PASS  $label_fmt  %s\n" "$1" "$2"; pass=$((pass+1)); }
fail_line() { printf "FAIL  $label_fmt  %s\n" "$1" "$2"; fail=$((fail+1)); }
skip_line() { printf "SKIP  $label_fmt  %s\n" "$1" "$2"; skip=$((skip+1)); }
detail()    { sed 's/^/          /'; }

# ── version ──────────────────────────────────────────────────────────────────

readme_v=$(sed -n '1s/^# SOD \([0-9][0-9.]*\).*/\1/p' README.md)
conf_v=$(sed -n 's/^release = "\([0-9][0-9.]*\)"/\1/p' docs/conf.py)
versions_v=$(sed -n 's/^## Version \([0-9][0-9.]*\).*/\1/p' VERSIONS.md | head -1)
banner_v=$(grep -ho "version [0-9][0-9]*\.[0-9][0-9]*" src/*.f90 | awk '{print $2}' | sort -u)

if [ -z "$readme_v" ] || [ -z "$conf_v" ] || [ -z "$versions_v" ]; then
    fail_line version "could not read a version from README.md, docs/conf.py or VERSIONS.md"
elif [ "$readme_v" != "$conf_v" ] || [ "$readme_v" != "$versions_v" ]; then
    fail_line version "README.md=$readme_v docs/conf.py=$conf_v VERSIONS.md=$versions_v"
elif [ "$(printf '%s\n' "$banner_v" | wc -l)" -ne 1 ] || [ "$banner_v" != "$readme_v" ]; then
    fail_line version "src/*.f90 startup banners disagree with $readme_v"
    printf '%s\n' "$banner_v" | detail
else
    pass_line version "$readme_v everywhere"
fi

# ── examples ─────────────────────────────────────────────────────────────────

missing=""
for d in examples/example[0-9][0-9]; do
    [ -d "$d" ] || continue
    name=$(basename "$d")
    grep -q "$name" README.md || missing="$missing $name"
done
if [ -n "$missing" ]; then
    fail_line examples "not described in README.md:$missing"
else
    pass_line examples "all $(ls -d examples/example[0-9][0-9] | wc -l | tr -d ' ') described in README.md"
fi

# ── scripts ──────────────────────────────────────────────────────────────────

# sod_common.sh is sourced, not called; the two maintenance scripts are not
# part of the documented workflow.
skip_scripts="sod_common.sh sod_run_tests.sh sod_check_docs.sh"
miss_readme=""; miss_docs=""
for f in bin/sod_*.sh; do
    s=$(basename "$f")
    case " $skip_scripts " in *" $s "*) continue ;; esac
    grep -q "$s" README.md            || miss_readme="$miss_readme $s"
    grep -rq "$s" docs/*.rst          || miss_docs="$miss_docs $s"
done
if [ -n "$miss_readme" ] || [ -n "$miss_docs" ]; then
    [ -n "$miss_readme" ] && fail_line scripts "missing from README.md:$miss_readme"
    [ -n "$miss_docs" ]   && fail_line scripts "missing from docs/:$miss_docs"
else
    pass_line scripts "every wrapper documented in both"
fi

# ── readme-toc ───────────────────────────────────────────────────────────────

toc_out=$(python3 - <<'PY'
import io, re, sys

text = io.open("README.md", encoding="utf-8").read()

def is_fence(line):
    s = line.strip()
    # a fence info string may not contain backticks, so ```cmd``` is inline code
    return s.startswith("```") and "`" not in s[3:]

def slug(title):
    s = re.sub(r"[^\w\s-]", "", title.lower(), flags=re.UNICODE)
    return s.replace(" ", "-")

fence = False
counts = {}
for line in text.split("\n"):
    if is_fence(line):
        fence = not fence
        continue
    if fence or not re.match(r"^#{1,6} ", line):
        continue
    a = slug(line.lstrip("#").strip())
    counts[a] = counts.get(a, 0) + 1

valid = set()
for a, n in counts.items():
    valid.add(a)
    valid.update(f"{a}-{i}" for i in range(1, n))

broken = sorted({m.group(1) for m in re.finditer(r"\]\(#([^)]+)\)", text)
                 if m.group(1) not in valid})
dups = sorted(a for a, n in counts.items() if n > 1)
links = len(re.findall(r"\]\(#", text))

if broken:
    print("BROKEN " + " ".join(broken))
elif dups:
    print("DUPS " + " ".join(dups))
else:
    print(f"OK {links}")
PY
) || toc_out="ERROR"
case "$toc_out" in
    OK*)     pass_line readme-toc "${toc_out#OK } in-page links resolve" ;;
    BROKEN*) fail_line readme-toc "unresolved anchors: ${toc_out#BROKEN }" ;;
    DUPS*)   fail_line readme-toc "duplicate heading slugs: ${toc_out#DUPS }" ;;
    *)       fail_line readme-toc "anchor check could not run" ;;
esac

# ── sphinx ───────────────────────────────────────────────────────────────────

sphinx=""
if command -v sphinx-build >/dev/null 2>&1; then
    sphinx="sphinx-build"
elif python3 -c "import sphinx" >/dev/null 2>&1; then
    sphinx="python3 -m sphinx"
fi
if [ -z "$sphinx" ]; then
    skip_line sphinx "sphinx-build not found (pip install -r docs/requirements.txt)"
else
    tmp=$(mktemp -d)
    if out=$($sphinx -b html -W --keep-going -q docs "$tmp" 2>&1); then
        pass_line sphinx "docs build clean with -W"
    else
        fail_line sphinx "docs build has warnings or errors"
        printf '%s\n' "$out" | head -20 | detail
    fi
    rm -rf "$tmp"
fi

# ── summary ──────────────────────────────────────────────────────────────────

echo
echo "passed: $pass   failed: $fail   skipped: $skip"
[ "$fail" -eq 0 ]
