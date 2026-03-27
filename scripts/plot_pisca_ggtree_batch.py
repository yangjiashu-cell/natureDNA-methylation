"""
Batch-call plot_pisca_ggtree.R on existing *_MCC.tree files under pisca_zenodo_fig4/ (or custom dir).

Skips if output PDF already exists unless --force. Does not fabricate trees — requires MCC from TreeAnnotator.

Requires: Rscript on PATH, R packages ape, ggtree, ggplot2.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
R_SCRIPT = ROOT / "scripts" / "plot_pisca_ggtree.R"
DEFAULT_DIR = ROOT / "figures" / "reproduced" / "pisca_zenodo_fig4"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, default=DEFAULT_DIR, help="Folder containing *_MCC.tree")
    ap.add_argument("--force", action="store_true", help="Overwrite existing PDFs")
    args = ap.parse_args()

    if not R_SCRIPT.exists():
        raise FileNotFoundError(R_SCRIPT)
    if shutil.which("Rscript") is None:
        print("Rscript not found on PATH. Install R and ensure Rscript is available.", file=sys.stderr)
        sys.exit(2)
    trees = sorted(args.dir.glob("*_MCC.tree"))
    if not trees:
        print(f"No *_MCC.tree under {args.dir}", file=sys.stderr)
        sys.exit(1)

    for t in trees:
        out = t.with_name(t.stem + "_ggtree.pdf")
        if out.exists() and not args.force:
            print("skip (exists):", out)
            continue
        cmd = ["Rscript", str(R_SCRIPT), str(t), str(out)]
        print("RUN:", " ".join(cmd))
        r = subprocess.run(cmd, cwd=str(ROOT))
        if r.returncode != 0:
            sys.exit(r.returncode)
    print("Done.")


if __name__ == "__main__":
    main()
