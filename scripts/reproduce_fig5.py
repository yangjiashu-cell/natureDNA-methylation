"""
Figure 5 (Nature paper layout) from MOESM9 (41586_2025_9374_MOESM9_ESM.xlsx).

  a — Univariate Cox: TTFT (blue) and OS (red) hazard ratios for five evolutionary
      variables (caption: EVOFLUx-inferred parameters).
  b — Kaplan–Meier: cumulative probability of first treatment (1 − S) by IGHV and
      high vs low inferred growth rate (theta), median split within each IGHV group.
  c — Multivariate Cox for TTFT: IGHV (U-CLL vs M-CLL), growth rate, TP53, age.

Outputs only PDFs: Fig5a.pdf, Fig5b.pdf, Fig5c.pdf (same names as paper panels a–c).

Methods: statsmodels PHReg (Breslow ties), Wald p-values; log-rank for KM comparisons
(approximate; see caption in original paper for Schoenfeld / Wald details).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.duration.hazard_regression import PHReg

ROOT = Path(__file__).resolve().parents[1]
XLSX = ROOT / "evoflux" / "41586_2025_9374_MOESM9_ESM.xlsx"
OUT = ROOT / "figures" / "reproduced"
OUT.mkdir(parents=True, exist_ok=True)

# Paper panel b: colour key (HTML + figure description)
COL_B = {
    "U-CLL High": "#8B1538",
    "U-CLL Low": "#F4A582",
    "M-CLL High": "#2E0854",
    "M-CLL Low": "#9B8FC8",
}


def _merge_5a_5c() -> pd.DataFrame:
    a = pd.read_excel(XLSX, sheet_name="Figure 5a")
    c = pd.read_excel(XLSX, sheet_name="Figure 5c")[
        ["PARTICIPANT_ID_ANONYMOUS", "DISEASE_SUBTYPE", "AGE_SAMPLING", "Genomics.Mutation_TP53"]
    ]
    m = a.merge(c, on="PARTICIPANT_ID_ANONYMOUS", how="inner")
    m = m[m["DISEASE_SUBTYPE"].notna()].copy()
    m["IGHV"] = m["DISEASE_SUBTYPE"].map({"mutated": "M-CLL", "unmutated": "U-CLL"})
    m = m[m["IGHV"].notna()]
    m["tp53_mut"] = (m["Genomics.Mutation_TP53"].astype(str) != "WT").astype(float)
    return m


def _cox_univariate(
    time_years: np.ndarray, event: np.ndarray, x: np.ndarray
) -> tuple[float, float, float, float] | None:
    """Return HR, CI_low, CI_high, Wald p-value; one continuous covariate."""
    x = np.asarray(x, dtype=float).ravel()
    ok = np.isfinite(time_years) & np.isfinite(event) & np.isfinite(x) & (time_years > 0)
    if ok.sum() < 40 or np.std(x[ok]) < 1e-12:
        return None
    t = time_years[ok]
    e = event[ok].astype(float)
    xv = x[ok].reshape(-1, 1)
    try:
        mod = PHReg(t, xv, status=e, ties="breslow")
        res = mod.fit(disp=False)
        beta = float(res.params[0])
        se = float(res.bse[0])
        hr = float(np.exp(beta))
        lo = float(np.exp(beta - 1.96 * se))
        hi = float(np.exp(beta + 1.96 * se))
        p = float(res.pvalues[0])
        return hr, lo, hi, p
    except Exception:
        return None


def _km_step(times_years: np.ndarray, events: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    t = np.asarray(times_years, dtype=float)
    e = np.asarray(events, dtype=int)
    order = np.argsort(t)
    t, e = t[order], e[order]
    uniq_event_times = np.sort(np.unique(t[e == 1]))
    xs = [0.0]
    ys = [1.0]
    S = 1.0
    for ti in uniq_event_times:
        n_risk = int(np.sum(t >= ti))
        d = int(np.sum((t == ti) & (e == 1)))
        if n_risk > 0 and d > 0:
            S *= 1.0 - d / n_risk
        xs.append(float(ti))
        ys.append(S)
    return np.array(xs), np.array(ys)


def _logrank_two_sample(t1: np.ndarray, e1: np.ndarray, t2: np.ndarray, e2: np.ndarray) -> float:
    """Two-sample log-rank p-value (chi-square df=1)."""
    t1, e1 = np.asarray(t1, float), np.asarray(e1, int)
    t2, e2 = np.asarray(t2, float), np.asarray(e2, int)
    times = np.concatenate([t1, t2])
    events = np.concatenate([e1, e2])
    groups = np.array([0] * len(t1) + [1] * len(t2), dtype=int)
    order = np.argsort(times)
    times, events, groups = times[order], events[order], groups[order]
    unique_times = np.unique(times[events == 1])
    oe = 0.0
    var = 0.0
    for ti in unique_times:
        at_risk_1 = int(np.sum((times >= ti) & (groups == 0)))
        at_risk_2 = int(np.sum((times >= ti) & (groups == 1)))
        at_risk = at_risk_1 + at_risk_2
        if at_risk < 2:
            continue
        d1 = int(np.sum((times == ti) & (events == 1) & (groups == 0)))
        d2 = int(np.sum((times == ti) & (events == 1) & (groups == 1)))
        d = d1 + d2
        if d == 0:
            continue
        e1_exp = at_risk_1 * d / at_risk
        oe += d1 - e1_exp
        var += (at_risk_1 * at_risk_2 * d * (at_risk - d)) / (at_risk**2 * (at_risk - 1)) if at_risk > 1 else 0.0
    if var <= 0:
        return 1.0
    chi2 = oe**2 / var
    return float(stats.chi2.sf(chi2, 1))


def _logrank_k_sample(groups_t: list[np.ndarray], groups_e: list[np.ndarray]) -> float:
    """K-sample log-rank p-value (chi-square df=k-1); groups_t/e same length lists."""
    k = len(groups_t)
    if k < 2:
        return 1.0
    all_t = np.concatenate([np.asarray(t, float) for t in groups_t])
    all_e = np.concatenate([np.asarray(e, int) for e in groups_e])
    groups = np.concatenate([np.full(len(groups_t[j]), j) for j in range(k)])
    order = np.argsort(all_t)
    all_t, all_e, groups = all_t[order], all_e[order], groups[order]
    uniq = np.unique(all_t[all_e == 1])
    O = np.zeros(k)
    E = np.zeros(k)
    cov = np.zeros((k, k))
    for ti in uniq:
        at_risk = np.array([np.sum((all_t >= ti) & (groups == j)) for j in range(k)], dtype=float)
        d = np.array([np.sum((all_t == ti) & (all_e == 1) & (groups == j)) for j in range(k)], dtype=float)
        dtot = float(d.sum())
        n = float(at_risk.sum())
        if dtot == 0 or n < 2:
            continue
        O += d
        E += at_risk * dtot / n
        fac = dtot * (n - dtot) / (n**2 * (n - 1.0))
        for i in range(k):
            for j in range(k):
                if i == j:
                    cov[i, j] += at_risk[i] * (n - at_risk[i]) * fac
                else:
                    cov[i, j] -= at_risk[i] * at_risk[j] * fac
    Z = O - E
    kk = k - 1
    Z1 = Z[:kk]
    V1 = cov[:kk, :kk]
    try:
        stat = float(Z1 @ np.linalg.solve(V1, Z1))
    except np.linalg.LinAlgError:
        return 1.0
    return float(stats.chi2.sf(stat, kk))


def _harrell_c_index(time: np.ndarray, event: np.ndarray, linpred: np.ndarray) -> float:
    """Concordance index for survival (risk = higher linpred = higher risk)."""
    t = np.asarray(time, float)
    e = np.asarray(event, int)
    r = np.asarray(linpred, float).ravel()
    n = len(t)
    num = 0
    den = 0
    for i in range(n):
        if e[i] != 1:
            continue
        for j in range(n):
            if i == j or t[j] < t[i]:
                continue
            if t[j] > t[i] or (t[j] == t[i] and e[j] == 0):
                den += 1
                if r[i] > r[j]:
                    num += 1
                elif r[i] == r[j]:
                    num += 0.5
    return float(num / den) if den > 0 else float("nan")


def fig5a() -> None:
    """Univariate survival analysis: TTFT (blue) and OS (red) hazard ratios."""
    df = _merge_5a_5c()
    vars_spec = [
        ("theta", "Growth rate"),
        ("Scancer", "Effective population size"),
        ("tau", "Patient's age at MRCA"),
        ("cancerAge", "Cancer age"),
        ("epiRate", "Mean epigenetic switching rate"),
    ]

    t_ttft = df["Clinics.TTFT_DAYS_SAMPLING"].to_numpy(dtype=float) / 365.25
    ev_ttft = df["Clinics.TTFT"].astype(float).to_numpy()
    t_os = df["Clinics.OS_DAYS_SAMPLING"].to_numpy(dtype=float) / 365.25
    ev_os = df["Clinics.OS"].astype(float).to_numpy()

    fig, ax = plt.subplots(figsize=(6.8, 5.2))
    n_rows = len(vars_spec)
    for i, (col, label) in enumerate(vars_spec):
        y0 = n_rows - 1 - i
        ax.axhspan(y0 - 0.48, y0 + 0.48, color="0.92" if i % 2 == 0 else "0.97", zorder=0)

        xraw = df[col].to_numpy(dtype=float)
        r_ttft = _cox_univariate(t_ttft, ev_ttft, xraw)
        r_os = _cox_univariate(t_os, ev_os, xraw)

        if r_ttft:
            hr, lo, hi, p = r_ttft
            ax.errorbar(
                hr,
                y0 + 0.12,
                xerr=[[hr - lo], [hi - hr]],
                fmt="s",
                color="#2c5aa0",
                ecolor="#2c5aa0",
                capsize=2.5,
                ms=5,
                zorder=3,
            )
            ax.text(
                hi * 1.1,
                y0 + 0.12,
                _format_p(p),
                fontsize=6.5,
                va="center",
                ha="left",
                color="#2c5aa0",
            )
        if r_os:
            hr, lo, hi, p = r_os
            ax.errorbar(
                hr,
                y0 - 0.12,
                xerr=[[hr - lo], [hi - hr]],
                fmt="o",
                color="#c0392b",
                ecolor="#c0392b",
                capsize=2.5,
                ms=5,
                zorder=3,
            )
            ax.text(
                hi * 1.1,
                y0 - 0.12,
                _format_p(p),
                fontsize=6.5,
                va="center",
                ha="left",
                color="#c0392b",
            )

        ax.text(-0.05, y0, label, ha="right", va="center", fontsize=8.5, transform=ax.get_yaxis_transform())

    ax.axvline(1.0, color="0.35", ls="--", lw=0.9, zorder=1)
    ax.set_xscale("log")
    ax.set_xlim(0.45, 5.5)
    ax.set_ylim(-0.6, n_rows - 0.4)
    ax.set_yticks([])
    ax.set_xlabel("Hazard ratio (95% CI)", fontsize=9)
    ax.set_title(
        "Fig. 5a | Univariate survival analysis of TTFT (blue) and OS (red)\n"
        "for evolutionary variables inferred via EVOFLUx (discovery cohort)",
        fontsize=9,
        pad=10,
    )

    # legend + sample counts
    ttft_ok = np.isfinite(t_ttft) & np.isfinite(ev_ttft)
    os_ok = np.isfinite(t_os) & np.isfinite(ev_os)
    n_ttft = int(ttft_ok.sum())
    ev_ttft_n = int(ev_ttft[ttft_ok].sum())
    n_os = int(os_ok.sum())
    ev_os_n = int(ev_os[os_ok].sum())
    ax.text(
        0.98,
        0.02,
        f"TTFT: n = {n_ttft}, events = {ev_ttft_n}\nOS: n = {n_os}, events = {ev_os_n}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.5,
        color="0.25",
    )
    from matplotlib.lines import Line2D

    leg = [
        Line2D([0], [0], color="#2c5aa0", marker="s", ls="None", label="TTFT"),
        Line2D([0], [0], color="#c0392b", marker="o", ls="None", label="OS"),
    ]
    ax.legend(handles=leg, loc="lower left", frameon=False, fontsize=8)
    sns.despine(ax=ax, left=True)
    fig.tight_layout()
    fig.savefig(OUT / "Fig5a.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _format_p(p: float) -> str:
    if not np.isfinite(p):
        return "—"
    if p < 1e-300:
        return r"$P \approx 0$"
    if p < 0.001:
        s = f"{p:.1e}"
        mant, exp = s.split("e")
        exp_i = int(exp)
        return rf"$P = {mant} \times 10^{{{exp_i}}}$"
    return rf"$P = {p:.4f}$"


def fig5b() -> None:
    """KM: cumulative probability of first treatment by IGHV × theta (median split)."""
    df = _merge_5a_5c().copy()
    df = df.dropna(subset=["theta", "Clinics.TTFT_DAYS_SAMPLING", "Clinics.TTFT", "DISEASE_SUBTYPE"])
    df = df[df["Clinics.TTFT_DAYS_SAMPLING"] > 0].copy()

    med_u = df[df["DISEASE_SUBTYPE"] == "unmutated"]["theta"].median()
    med_m = df[df["DISEASE_SUBTYPE"] == "mutated"]["theta"].median()

    def strat(row: pd.Series) -> str:
        if row["DISEASE_SUBTYPE"] == "unmutated":
            h = row["theta"] >= med_u
            return "U-CLL High" if h else "U-CLL Low"
        h = row["theta"] >= med_m
        return "M-CLL High" if h else "M-CLL Low"

    df["stratum"] = df.apply(strat, axis=1)

    t = df["Clinics.TTFT_DAYS_SAMPLING"].to_numpy(dtype=float) / 365.25
    ev = df["Clinics.TTFT"].astype(int).to_numpy()

    fig, ax = plt.subplots(figsize=(5.5, 6.2))
    order = ["U-CLL High", "U-CLL Low", "M-CLL High", "M-CLL Low"]
    for lab in order:
        m = df["stratum"] == lab
        if not m.any():
            continue
        tx, sx = _km_step(t[m.to_numpy()], ev[m.to_numpy()])
        # cumulative probability of treatment = 1 - S(t) (still treatment-free)
        cum = 1.0 - sx
        ax.step(tx, cum, where="post", color=COL_B[lab], lw=1.8, label=f"{lab} (n = {m.sum()})")

    ax.set_xlim(0, 12)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Time from sampling (years)", fontsize=9)
    ax.set_ylabel("Probability of treatment", fontsize=9)
    ax.set_title(
        "Fig. 5b | Kaplan–Meier curves for TTFT (high vs low growth rate within IGHV)\n"
        "Cumulative probability of first treatment; median split on θ",
        fontsize=9,
        pad=8,
    )
    ax.legend(loc="upper left", frameon=False, fontsize=7.5)

    # Log-rank: U-CLL high vs low, M-CLL high vs low
    u = df["DISEASE_SUBTYPE"] == "unmutated"
    m = df["DISEASE_SUBTYPE"] == "mutated"
    p_u = _logrank_two_sample(
        t[u & (df["stratum"] == "U-CLL High")],
        ev[u & (df["stratum"] == "U-CLL High")],
        t[u & (df["stratum"] == "U-CLL Low")],
        ev[u & (df["stratum"] == "U-CLL Low")],
    )
    p_m = _logrank_two_sample(
        t[m & (df["stratum"] == "M-CLL High")],
        ev[m & (df["stratum"] == "M-CLL High")],
        t[m & (df["stratum"] == "M-CLL Low")],
        ev[m & (df["stratum"] == "M-CLL Low")],
    )
    groups_t: list[np.ndarray] = []
    groups_e: list[np.ndarray] = []
    for lab in order:
        mask = (df["stratum"] == lab).to_numpy()
        if not mask.any():
            continue
        groups_t.append(t[mask])
        groups_e.append(ev[mask])
    p_all = _logrank_k_sample(groups_t, groups_e)
    ax.text(0.02, 0.98, f"Four groups, log-rank: {_format_p(p_all)}", transform=ax.transAxes, fontsize=7, va="top")
    ax.text(0.02, 0.93, f"U-CLL high vs low: {_format_p(p_u)}", transform=ax.transAxes, fontsize=7, va="top")
    ax.text(0.02, 0.88, f"M-CLL high vs low: {_format_p(p_m)}", transform=ax.transAxes, fontsize=7, va="top")
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(OUT / "Fig5b.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fig5c() -> None:
    """Multivariate Cox for TTFT: IGHV, theta, TP53, age."""
    df = _merge_5a_5c()
    df = df.dropna(
        subset=[
            "Clinics.TTFT_DAYS_SAMPLING",
            "Clinics.TTFT",
            "theta",
            "AGE_SAMPLING",
            "DISEASE_SUBTYPE",
        ]
    )
    df = df[df["Clinics.TTFT_DAYS_SAMPLING"] > 0].copy()
    df["IGHV_U"] = (df["DISEASE_SUBTYPE"] == "unmutated").astype(float).values  # U-CLL vs M-CLL ref

    t = df["Clinics.TTFT_DAYS_SAMPLING"].to_numpy(dtype=float) / 365.25
    ev = df["Clinics.TTFT"].astype(float).to_numpy()
    X = np.column_stack(
        [
            df["IGHV_U"].to_numpy(float),
            df["theta"].to_numpy(float),
            df["tp53_mut"].to_numpy(float),
            df["AGE_SAMPLING"].to_numpy(float),
        ]
    )
    labels = ["IGHV U-CLL", "Growth rate", "TP53", "Age"]

    mod = PHReg(t, X, status=ev, ties="breslow")
    res = mod.fit(disp=False)
    beta = np.asarray(res.params, dtype=float).ravel()
    hrs = np.exp(beta)
    ci_beta = np.asarray(res.conf_int())
    ci_lo, ci_hi = np.exp(ci_beta[:, 0]), np.exp(ci_beta[:, 1])
    ps = res.pvalues

    # Global Wald test: H0: all coefficients zero
    try:
        global_p = float(res.wald_test(np.eye(len(beta)), scalar=True).pvalue)
    except Exception:
        global_p = float("nan")

    # C-index
    linpred = X @ beta
    cidx = _harrell_c_index(t, ev.astype(int), linpred)

    n = len(df)
    ev_n = int(ev.sum())

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    y = np.arange(len(labels))
    ax.axvline(1.0, color="0.35", ls="--", lw=0.9)
    for i in range(len(labels)):
        ax.axhspan(i - 0.45, i + 0.45, color="0.92" if i % 2 == 0 else "0.97", zorder=0)
    ax.errorbar(
        hrs,
        y,
        xerr=[hrs - ci_lo, ci_hi - hrs],
        fmt="s",
        color="#2c5aa0",
        ecolor="#2c5aa0",
        capsize=2.5,
        ms=5,
        zorder=3,
    )
    for i, lab in enumerate(labels):
        ax.text(
            ci_hi[i] * 1.08,
            i,
            _format_p(float(np.asarray(ps).ravel()[i])),
            fontsize=6.5,
            va="center",
            ha="left",
            color="#2c5aa0",
        )

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xscale("log")
    ax.set_xlim(0.45, 8)
    ax.set_xlabel("Hazard ratio (95% CI) for TTFT", fontsize=9)
    ax.set_title(
        "Fig. 5c | Multivariate Cox regression model of TTFT\n"
        "(IGHV U-CLL vs M-CLL, growth rate, TP53, age at sampling)",
        fontsize=9,
        pad=10,
    )
    ax.text(
        0.98,
        0.02,
        f"n = {n}\nEvents = {ev_n}\nC-index ≈ {cidx:.2f}\nGlobal {_format_p(global_p)}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.5,
        color="0.25",
    )
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(OUT / "Fig5c.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "pdf.fonttype": 42,
        }
    )
    fig5a()
    print(f"[Fig5a] {OUT / 'Fig5a.pdf'}")
    fig5b()
    print(f"[Fig5b] {OUT / 'Fig5b.pdf'}")
    fig5c()
    print(f"[Fig5c] {OUT / 'Fig5c.pdf'}")


if __name__ == "__main__":
    main()
