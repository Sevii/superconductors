#!/usr/bin/env bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
ROOT="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$ROOT/certificates" "$ROOT/data"
find "$ROOT" -type d -name __pycache__ -prune -exec rm -rf {} +
find "$ROOT" -type f -name '*.pyc' -delete

cd "$ROOT/scripts"
python3 verify_hminusone_transfer.py > "$ROOT/certificates/HMINUSONE_TRANSFER_CERTIFICATE.txt"
python3 verify_two_control_overlap.py > "$ROOT/certificates/TWO_CONTROL_OVERLAP_CERTIFICATE.txt"
python3 audit_order8_all_fillings.py --json "$ROOT/data/order8_all_fillings.json" > "$ROOT/certificates/ORDER8_ALL_FILLINGS_OUTPUT.txt"
python3 analyze_order8_tuning.py > "$ROOT/certificates/ORDER8_TUNING_OUTPUT.txt"
python3 analyze_kramers_order8.py > "$ROOT/certificates/KRAMERS_ORDER8_OUTPUT.txt"
python3 verify_order8_structure.py > "$ROOT/certificates/ORDER8_STRUCTURE_CERTIFICATE.txt"
python3 verify_kramers_feshbach.py > "$ROOT/certificates/KRAMERS_FESHBACH_CERTIFICATE.txt"
python3 verify_dark_agp_branch.py > "$ROOT/certificates/DARK_AGP_BRANCH_CERTIFICATE.txt"

[ -f order8_tuning_analysis.json ] && mv -f order8_tuning_analysis.json "$ROOT/data/"
[ -f kramers_order8_analysis.json ] && mv -f kramers_order8_analysis.json "$ROOT/data/"

cd "$ROOT/paper"
rm -f microscopic_transfer_exact_meissner.{aux,fdb_latexmk,fls,log,out,toc}
if command -v latexmk >/dev/null 2>&1; then
  latexmk -pdf -interaction=nonstopmode -halt-on-error microscopic_transfer_exact_meissner.tex > BUILD_LATEXMK.log
  # Retain a clean final pass rather than latexmk's expected first-pass warnings.
  pdflatex -interaction=nonstopmode -halt-on-error microscopic_transfer_exact_meissner.tex > BUILD_LOG.txt
  rm -f BUILD_LATEXMK.log
else
  pdflatex -interaction=nonstopmode -halt-on-error microscopic_transfer_exact_meissner.tex > BUILD_PASS_1.log
  pdflatex -interaction=nonstopmode -halt-on-error microscopic_transfer_exact_meissner.tex > BUILD_PASS_2.log
  pdflatex -interaction=nonstopmode -halt-on-error microscopic_transfer_exact_meissner.tex > BUILD_LOG.txt
  rm -f BUILD_PASS_1.log BUILD_PASS_2.log
fi
rm -f microscopic_transfer_exact_meissner.{aux,fdb_latexmk,fls,log,out,toc}

cd "$ROOT"
find "$ROOT" -type d -name __pycache__ -prune -exec rm -rf {} +
find "$ROOT" -type f -name '*.pyc' -delete
python3 scripts/verify_release_hygiene.py > certificates/RELEASE_HYGIENE_CERTIFICATE.txt

printf '\nCertificate summary\n'
for f in "$ROOT"/certificates/*.txt; do
  printf '%-42s ' "$(basename "$f")"
  tail -n 1 "$f"
done
printf '%-42s ' "PDF pages"
pdfinfo "$ROOT/paper/microscopic_transfer_exact_meissner.pdf" | awk '/^Pages:/ {print $2}'
