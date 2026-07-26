#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

python check_multiplet_uv_completion.py \
  > multiplet_uv_completion_certificate.txt
python check_electron_only_parity_channels.py \
  > electron_only_parity_channels_certificate.txt
python check_all_filling_active_remote_completion.py \
  > all_filling_active_remote_completion_certificate.txt
python check_global_lattice_downfolding.py \
  > global_lattice_downfolding_certificate.txt
python check_kramers_closed_shell_embedding.py \
  > kramers_closed_shell_embedding_certificate.txt

for certificate in *_certificate.txt; do
  if ! grep -qE 'OVERALL: PASS|All microscopic multiplet checks passed' "$certificate"; then
    echo "Certificate did not report PASS: $certificate" >&2
    exit 1
  fi
done

echo "All five certificates: PASS"
