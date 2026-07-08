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
