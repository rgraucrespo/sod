#!/bin/bash
# Shared helpers for SOD wrapper scripts. Source this file; do not execute it.

sod_find_project_root() {
  local dir
  dir="$(cd "${1:-$PWD}" && pwd -P)" || return 1

  while :; do
    if [ -f "$dir/INSOD" ] && [ -f "$dir/SGO" ]; then
      printf '%s\n' "$dir"
      return 0
    fi
    [ "$dir" = "/" ] && return 1
    dir="$(dirname "$dir")"
  done
}

sod_require_project_root() {
  local root
  if root="$(sod_find_project_root "${1:-$PWD}")"; then
    printf '%s\n' "$root"
    return 0
  fi
  echo "Error: could not locate SODPROJECT/ from $(pwd)." >&2
  echo "       Searched upward for a directory containing INSOD and SGO." >&2
  return 1
}

sod_read_insod_filer() {
  local insod_path last_data_line
  insod_path="${1:-INSOD}"

  if [ ! -f "$insod_path" ]; then
    echo "Error: INSOD file not found: $insod_path" >&2
    return 1
  fi

  last_data_line="$(
    awk '
      {
        line = $0
        sub(/^[[:space:]]+/, "", line)
        sub(/[[:space:]]+$/, "", line)
        if (line != "" && substr(line, 1, 1) != "#") last = line
      }
      END {
        if (last == "") exit 1
        print last
      }
    ' "$insod_path"
  )" || {
    echo "Error: could not locate FILER line in $insod_path" >&2
    return 1
  }

  if [[ ! "$last_data_line" =~ ^-?[0-9]+$ ]]; then
    echo "Error: could not parse FILER from INSOD tail line: $last_data_line" >&2
    return 1
  fi

  printf '%s\n' "$last_data_line"
}

sod_parse_model_args() {
  SOD_MODEL_ARGS=()

  while [[ $# -gt 0 ]]; do
    case "$1" in
      -model)
        if [[ $# -lt 2 || -z "$2" ]]; then
          echo "Error: -model requires a filename argument." >&2
          return 1
        fi
        SOD_MODEL_ARGS+=("-model" "$2")
        shift 2
        ;;
      *)
        echo "Error: unexpected argument: $1" >&2
        return 1
        ;;
    esac
  done
}

sod_find_enclosing_level_name() {
  local root dir base
  root="$(cd "$1" && pwd -P)" || return 1
  dir="$(cd "${2:-$PWD}" && pwd -P)" || return 1

  while [ "$dir" != "$root" ] && [ "$dir" != "/" ]; do
    base="$(basename "$dir")"
    if [[ "$base" =~ ^n[0-9]+(_[0-9]+)*$ ]]; then
      printf '%s\n' "$base"
      return 0
    fi
    dir="$(dirname "$dir")"
  done
  return 1
}

sod_find_ancestor_with_file() {
  local root dir filename
  root="$(cd "$1" && pwd -P)" || return 1
  dir="$(cd "${2:-$PWD}" && pwd -P)" || return 1
  filename="$3"

  while [ "$dir" != "/" ]; do
    if [ -f "$dir/$filename" ]; then
      printf '%s\n' "$dir"
      return 0
    fi
    [ "$dir" = "$root" ] && return 1
    dir="$(dirname "$dir")"
  done
  return 1
}

sod_level_dir_by_number() {
  local root level d name digits
  root="$1"
  level="$2"

  for d in "$root"/n*/; do
    [ -d "$d" ] || continue
    name="$(basename "${d%/}")"
    if [[ "$name" =~ ^n0*([0-9]+)$ ]]; then
      digits="${BASH_REMATCH[1]}"
      if [ "$((10#$digits))" -eq "$level" ]; then
        printf '%s\n' "$name"
        return 0
      fi
    fi
  done
  return 1
}

# Read the configuration count from the first line of a version 3 ENSEMBLE.
#
# Exit status: 0 with the count on stdout; 2 if the file is missing or is not
# a version 3 ENSEMBLE. Version 2 files are no longer supported; regenerate
# them with sod_comb.sh or sod_random.sh.
sod_ensemble_config_count() {
  local ensemble_path count
  ensemble_path="$1"

  if [ ! -f "$ensemble_path" ]; then
    echo "Error: ENSEMBLE file not found: $ensemble_path" >&2
    return 2
  fi

  # Line 1 of a v3 ENSEMBLE is "<type> ensemble[ (T K)]: N configurations; ...",
  # so the count is the integer token immediately before "configurations".
  count="$(
    awk '
      /^[[:space:]]*$/ { next }
      {
        if ($1 !~ /^#/) {
          for (i = 2; i <= NF; i++) {
            if ($i ~ /^configurations[;:]?$/ && $(i-1) ~ /^[0-9]+$/) {
              print $(i-1)
              break
            }
          }
        }
        exit
      }
    ' "$ensemble_path"
  )"
  if [[ ! "$count" =~ ^[0-9]+$ ]] || [ "$count" -lt 1 ]; then
    echo "Error: $ensemble_path is not a version 3 ENSEMBLE." >&2
    echo "       Expected a first line such as" >&2
    echo "       'Enumerated ensemble: 71 configurations; sum_degeneracies = 4096'." >&2
    echo "       Version 2 ENSEMBLE files are no longer supported; regenerate the" >&2
    echo "       level with sod_comb.sh or sod_random.sh." >&2
    return 2
  fi
  printf '%s\n' "$count"
}

# Does level_dir hold output_name for every configuration 1..N of its ENSEMBLE?
#
# Exit status: 0 complete; 1 incomplete (missing, duplicated or out-of-range
# cYY indices); 2 the ENSEMBLE could not be read. Sets
# SOD_ENSEMBLE_CONFIG_COUNT and SOD_COMPLETE_OUTPUT_COUNT for statuses 0 and 1.
sod_level_outputs_complete() {
  local level_dir output_name expected dir name digits index available invalid
  local -A present=()
  level_dir="$1"
  output_name="$2"

  expected="$(sod_ensemble_config_count "$level_dir/ENSEMBLE")" || return 2
  available=0
  invalid=0
  for dir in "$level_dir"/c[0-9]*/; do
    [ -d "$dir" ] || continue
    name="$(basename "${dir%/}")"
    digits="${name#c}"
    # Only all-digit cYY names are configuration directories, matching
    # sodpaths.config_dirs(); c001_backup and the like are not.
    [[ "$digits" =~ ^[0-9]+$ ]] || continue
    [ -f "$dir/$output_name" ] || continue
    index="$((10#$digits))"
    if [ "$index" -lt 1 ] || [ "$index" -gt "$expected" ] || [ -n "${present[$index]+x}" ]; then
      invalid=1
    else
      present[$index]=1
      available=$((available+1))
    fi
  done

  SOD_ENSEMBLE_CONFIG_COUNT="$expected"
  SOD_COMPLETE_OUTPUT_COUNT="$available"
  [ "$available" -eq "$expected" ] && [ "$invalid" -eq 0 ]
}

# Print the cYY configuration directory names of a level directory (default: the
# current directory), one per line, in numeric index order. Only all-digit cYY
# names qualify, matching sodpaths.config_dirs(); c001_backup, cSGO and the like
# are skipped. Names are printed bare, without the leading directory, so callers
# compose the path they need.
sod_config_dirs_in_order() {
  local base dir name digits
  base="${1:-.}"
  for dir in "${base%/}"/c[0-9]*/; do
    [ -d "$dir" ] || continue
    name="${dir%/}"
    name="${name##*/}"
    digits="${name#c}"
    [[ "$digits" =~ ^[0-9]+$ ]] || continue
    printf '%d\t%s\n' "$((10#$digits))" "$name"
  done | sort -n -k1,1 | cut -f2
}
