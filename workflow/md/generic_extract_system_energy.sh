#!/usr/bin/env bash
set -euo pipefail

GMX_BIN="${GMX_BIN:-/usr/local/gmx2025.4/bin/gmx}"

Edr=""
OutXvg=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --edr) Edr="$2"; shift 2 ;;
        --out-xvg) OutXvg="$2"; shift 2 ;;
        *)
            printf 'Unknown argument: %s\n' "$1" >&2
            exit 1
            ;;
    esac
done

[[ -f "$Edr" ]] || { printf 'Missing EDR: %s\n' "$Edr" >&2; exit 1; }
[[ -n "$OutXvg" ]] || { printf '--out-xvg is required\n' >&2; exit 1; }

mkdir -p "$(dirname "$OutXvg")"
printf 'Potential\n0\n' | "$GMX_BIN" energy -f "$Edr" -o "$OutXvg" >/dev/null
