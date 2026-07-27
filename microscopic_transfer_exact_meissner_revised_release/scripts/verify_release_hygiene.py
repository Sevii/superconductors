#!/usr/bin/env python3
"""Deterministic release-hygiene certificate.

Checks all shipped text-bearing files for embedded C0 control characters other
than LF/CR and rejects Python bytecode/cache artifacts. This specifically guards
against escaped ``\\t`` and ``\\f`` sequences becoming literal TAB/form-feed
characters in generated TeX or Markdown.
"""
from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
TEXT_SUFFIXES = {".tex", ".md", ".py", ".sh", ".txt", ".json"}
ALLOWED_CONTROLS = {10, 13}  # LF, CR


def main() -> int:
    bad_controls: list[tuple[str, int, int]] = []
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        for offset, byte in enumerate(path.read_bytes()):
            if byte < 32 and byte not in ALLOWED_CONTROLS:
                bad_controls.append((str(path.relative_to(ROOT)), offset, byte))

    bytecode = sorted(
        str(path.relative_to(ROOT))
        for path in ROOT.rglob("*")
        if path.is_file() and path.suffix == ".pyc"
    )
    caches = sorted(
        str(path.relative_to(ROOT))
        for path in ROOT.rglob("__pycache__")
        if path.is_dir()
    )

    print("RELEASE HYGIENE CERTIFICATE")
    print("=" * 72)
    print(f"text files scanned: {sum(1 for p in ROOT.rglob('*') if p.is_file() and p.suffix.lower() in TEXT_SUFFIXES)}")
    print(f"embedded disallowed C0 controls: {len(bad_controls)}")
    for rel, offset, byte in bad_controls[:20]:
        print(f"  {rel}: byte offset {offset}, value 0x{byte:02x}")
    print(f"Python .pyc files: {len(bytecode)}")
    for rel in bytecode[:20]:
        print(f"  {rel}")
    print(f"__pycache__ directories: {len(caches)}")
    for rel in caches[:20]:
        print(f"  {rel}")

    passed = not bad_controls and not bytecode and not caches
    print("OVERALL:", "PASS" if passed else "FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
