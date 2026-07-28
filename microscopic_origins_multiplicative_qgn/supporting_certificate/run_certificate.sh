#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
python3 "$HERE/check_unscreened_agp_leakage.py" --output-dir "$HERE"
