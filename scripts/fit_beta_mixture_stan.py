"""
Optional Stan (CmdStan) fit for per-sample 3-state methylation coding.

Requires: cmdstanpy + CmdStan installed; set CMDSTAN environment variable.

This does NOT fabricate results: if CmdStan is missing, the script exits with an error.

Usage:
  python scripts/fit_beta_mixture_stan.py --sample SCLL-545 \\
    --beta-csv data/beta_fcpgs.csv --out figures/reproduced/pisca_zenodo_fig4/SCLL-545_stan_states.txt

Batch (Fig.4 Zenodo cases) -> CSV for --discretize external:
  python scripts/fit_beta_mixture_stan.py --case case12 --out-csv figures/reproduced/pisca_zenodo_fig4/case12_stan_states_matrix.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
STAN = ROOT / "stan" / "beta_mixture_3state.stan"

CASE12_SAMPLES = ["SCLL-545", "SCLL-546", "SCLL-547", "SCLL-548"]
CASE19_SAMPLES = ["SCLL-493", "SCLL-494", "SCLL-531", "SCLL-532", "SCLL-533"]


def _fit_one_sample(
    beta: pd.DataFrame,
    sample: str,
    *,
    chains: int,
    iter_warmup: int,
    iter_sampling: int,
) -> list[int]:
    try:
        import cmdstanpy
    except ImportError as e:
        raise SystemExit("Install cmdstanpy and CmdStan to use Stan: " + str(e)) from e

    if not STAN.exists():
        raise FileNotFoundError(STAN)
    if sample not in beta.columns:
        raise SystemExit(f"Sample {sample!r} not in beta matrix")
    y = np.asarray(beta[sample], dtype=float)
    y = np.clip(y, 1e-4, 1.0 - 1e-4)

    sm = cmdstanpy.CmdStanModel(stan_file=str(STAN))
    fit = sm.sample(
        data={"N": int(y.size), "y": y.tolist()},
        chains=chains,
        iter_warmup=iter_warmup,
        iter_sampling=iter_sampling,
        show_progress=True,
        show_console=False,
    )
    from scipy.stats import beta as beta_dist

    draws = fit.draws_pd()
    means: list[float] = []
    for k in (1, 2, 3):
        a = float(draws[f"alpha[{k}]"].median())
        b = float(draws[f"beta[{k}]"].median())
        means.append(a / (a + b))
    order = np.argsort(means)
    inv = np.empty(3, dtype=int)
    inv[order] = np.arange(3)

    ta = [float(draws[f"alpha[{k}]"].median()) for k in (1, 2, 3)]
    tb = [float(draws[f"beta[{k}]"].median()) for k in (1, 2, 3)]
    tt = [float(draws[f"theta[{k}]"].median()) for k in (1, 2, 3)]
    s = sum(tt)
    tt = [t / s for t in tt]

    states: list[int] = []
    for val in y:
        logp = [np.log(tt[k]) + beta_dist.logpdf(val, ta[k], tb[k]) for k in range(3)]
        k_hat = int(np.argmax(logp))
        states.append(int(inv[k_hat]))
    return states


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", default=None, help="Column name in beta_fcpgs.csv")
    ap.add_argument(
        "--case",
        choices=("case12", "case19"),
        default=None,
        help="Fit all longitudinal samples for Fig.4 case (writes matrix if --out-csv set).",
    )
    ap.add_argument("--beta-csv", type=Path, default=ROOT / "data" / "beta_fcpgs.csv")
    ap.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output text: one digit 0/1/2 per CpG row order (single --sample)",
    )
    ap.add_argument(
        "--out-csv",
        type=Path,
        default=None,
        help="Samples × CpG matrix (integers 0–2) for reproduce_pisca_zenodo_fig4.py --discretize external",
    )
    ap.add_argument("--chains", type=int, default=2)
    ap.add_argument("--iter-warmup", type=int, default=500)
    ap.add_argument("--iter-sampling", type=int, default=500)
    args = ap.parse_args()

    if args.case and args.sample:
        raise SystemExit("Use either --case or --sample, not both")
    if not args.case and not args.sample:
        raise SystemExit("Provide --sample or --case")
    if args.sample and not args.out:
        raise SystemExit("--sample requires --out")
    if args.case and not args.out_csv:
        raise SystemExit("--case requires --out-csv")

    beta = pd.read_csv(args.beta_csv, index_col=0)

    if args.case:
        samples = CASE12_SAMPLES if args.case == "case12" else CASE19_SAMPLES
        rows: dict[str, list[int]] = {}
        for sid in samples:
            print(f"Stan fit: {sid} ...")
            rows[sid] = _fit_one_sample(
                beta,
                sid,
                chains=args.chains,
                iter_warmup=args.iter_warmup,
                iter_sampling=args.iter_sampling,
            )
        cpg_ids = beta.index.astype(str)
        mat = pd.DataFrame([rows[s] for s in samples], index=samples, columns=cpg_ids)
        args.out_csv.parent.mkdir(parents=True, exist_ok=True)
        mat.to_csv(args.out_csv)
        print(f"Wrote {args.out_csv} shape {mat.shape}")
        return

    states = _fit_one_sample(
        beta,
        args.sample,
        chains=args.chains,
        iter_warmup=args.iter_warmup,
        iter_sampling=args.iter_sampling,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("".join(str(s) for s in states), encoding="utf-8")
    print(f"Wrote {args.out} ({len(states)} sites)")


if __name__ == "__main__":
    main()
