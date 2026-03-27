"""
Run BEAST 1.10.x on a prepared XML (SNV / binary pipeline).

PISCA is documented as incompatible with BEAST v1.10; SNV phylogeny uses standard BEAST models.

Does not fabricate results: exits non-zero if BEAST 1.10 is missing or XML is missing.

Usage:
  python scripts/run_snv_beast10.py --xml figures/reproduced/snv_beast10_prep/case_1.xml
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _find_beast10() -> Path | None:
    candidates = [
        ROOT / "tools" / "BEAST v1.10.0",
        ROOT / "tools" / "BEAST v1.10.4",
        ROOT / "tools" / "BEAST v1.10.5",
        Path(r"D:\tools\BEASTv1.10.5"),
        Path(r"C:\tools\BEASTv1.10.5"),
    ]
    for c in candidates:
        cmd = c / "bin" / "beast.cmd"
        if cmd.exists() and (c / "lib" / "beast.jar").exists():
            return cmd
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xml", type=Path, required=True)
    ap.add_argument(
        "--overwrite",
        action="store_true",
        help="Pass -overwrite to BEAST (recommended for scripted re-runs).",
    )
    args = ap.parse_args()

    beast = _find_beast10()
    if beast is None:
        print("BEAST 1.10.x not found under tools/ or D:\\tools\\BEASTv1.10.5", file=sys.stderr)
        sys.exit(2)
    if not args.xml.exists():
        print(f"XML not found: {args.xml}", file=sys.stderr)
        sys.exit(2)

    cmd = [str(beast)]
    if args.overwrite:
        cmd.append("-overwrite")
    cmd.append(str(args.xml))
    print("Running:", " ".join(cmd))
    r = subprocess.run(cmd, cwd=str(args.xml.parent))
    sys.exit(r.returncode)


if __name__ == "__main__":
    main()
