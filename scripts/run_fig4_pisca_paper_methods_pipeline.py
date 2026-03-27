"""
Orchestrate Fig.4 PISCA steps aligned with Nature 2025 Methods (Zenodo β + ages):

  1) Stan 3-state Beta mixture per sample → CSV matrices (scripts/fit_beta_mixture_stan.py)
  2) BEAST 1.8.4 + PISCA XML (scripts/reproduce_pisca_zenodo_fig4.py --discretize external)
  3) Optional: run BEAST / TreeAnnotator (same CLI flags as reproduce script)

Author repo CalumGabbutt/evoflux ships EVOFLUx inference only; it does not include PISCA/BEAST/Stan
phylogeny code. This pipeline follows the paper Methods + PISCAv1.1 template already in tools/.

Does not fabricate results: if CmdStan is missing, step (1) fails unless --skip-stan and existing CSVs.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIT_STAN = ROOT / "scripts" / "fit_beta_mixture_stan.py"
REPRO = ROOT / "scripts" / "reproduce_pisca_zenodo_fig4.py"
OUT = ROOT / "figures" / "reproduced" / "pisca_zenodo_fig4"

CASE_MAP = {
    "case12": OUT / "case12_stan_states_matrix.csv",
    "case19": OUT / "case19_stan_states_matrix.csv",
}


def _run(cmd: list[str]) -> int:
    print("RUN:", " ".join(cmd))
    return subprocess.call(cmd, cwd=str(ROOT))


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Stan → external CSV → PISCA XML for Fig.4 (paper Methods)."
    )
    ap.add_argument("--cases", default="case12,case19", help="Comma-separated: case12, case19")
    ap.add_argument(
        "--skip-stan",
        action="store_true",
        help="Skip CmdStan step; CSVs must already exist under figures/reproduced/pisca_zenodo_fig4/",
    )
    ap.add_argument(
        "--paper-chain",
        action="store_true",
        help="Pass through to reproduce script (1e8 generations; very long).",
    )
    ap.add_argument(
        "--also-write-exponential-xml",
        action="store_true",
        help="Pass through: write <stem>_exponential.xml for path sampling / marginal likelihood.",
    )
    ap.add_argument(
        "--write-path-sampling-appendix",
        action="store_true",
        help="Write path_sampling_beast18_notes.txt (no MCMC).",
    )
    ap.add_argument("--skip-beast", action="store_true", help="XML + CSV only (no BEAST).")
    ap.add_argument("--no-annotator", action="store_true")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)

    keys = [k.strip() for k in args.cases.split(",") if k.strip()]
    for key in keys:
        if key not in CASE_MAP:
            raise SystemExit(f"Unknown case {key!r}; choose from {list(CASE_MAP)}")

    if not args.skip_stan:
        for key in keys:
            csv_path = CASE_MAP[key]
            code = _run(
                [
                    sys.executable,
                    str(FIT_STAN),
                    "--case",
                    key,
                    "--out-csv",
                    str(csv_path),
                ]
            )
            if code != 0:
                raise SystemExit(f"fit_beta_mixture_stan failed for {key} (exit {code})")

    for key in keys:
        csv_path = CASE_MAP[key]
        if not csv_path.exists():
            raise FileNotFoundError(
                f"Missing {csv_path}. Run without --skip-stan or place Stan output CSV here."
            )
        cmd = [
            sys.executable,
            str(REPRO),
            "--cases",
            key,
            "--discretize",
            "external",
            "--external-states-csv",
            str(csv_path),
        ]
        if args.paper_chain:
            cmd.append("--paper-chain")
        if args.also_write_exponential_xml:
            cmd.append("--also-write-exponential-xml")
        if args.write_path_sampling_appendix:
            cmd.append("--write-path-sampling-appendix")
        if args.skip_beast:
            cmd.append("--skip-beast")
        if args.no_annotator:
            cmd.append("--no-annotator")
        code = _run(cmd)
        if code != 0:
            raise SystemExit(f"reproduce_pisca_zenodo_fig4 failed for {key} (exit {code})")

    print("Pipeline finished. Outputs under", OUT)


if __name__ == "__main__":
    main()
