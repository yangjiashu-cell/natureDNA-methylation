from __future__ import annotations

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BEAST_BIN = ROOT / "tools" / "BEAST v1.8.4" / "bin"
TREEANNOTATOR = BEAST_BIN / "treeannotator.cmd"
OUT_DIR = ROOT / "figures" / "reproduced" / "pisca_moesm8_fallback"


def run_mcc(
    stem: str,
    *,
    burnin_trees: int = 5,
    heights: str = "mean",
) -> Path:
    """Run TreeAnnotator (BEAST 1.8) to produce an annotated MCC tree from posterior .trees."""
    trees = OUT_DIR / f"{stem}.trees"
    out = OUT_DIR / f"{stem}_MCC.tree"
    if not trees.exists():
        raise FileNotFoundError(trees)
    if not TREEANNOTATOR.exists():
        raise FileNotFoundError(TREEANNOTATOR)

    cmd = [
        str(TREEANNOTATOR),
        "-burninTrees",
        str(burnin_trees),
        "-heights",
        heights,
        str(trees),
        str(out),
    ]
    completed = subprocess.run(cmd, cwd=str(OUT_DIR), check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"TreeAnnotator failed with exit code {completed.returncode}")
    return out


def main() -> None:
    for stem in ("fig4a_case12_moesm8_fallback", "fig4b_case19_moesm8_fallback"):
        p = run_mcc(stem)
        print(f"WROTE {p}")


if __name__ == "__main__":
    main()
