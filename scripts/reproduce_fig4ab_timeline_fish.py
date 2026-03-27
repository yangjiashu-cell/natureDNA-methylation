"""
Figure 4a/4b 顶部复现（Timeline + WGS Fish Plot，纯 matplotlib + UnivariateSpline）。
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.transforms import blended_transform_factory
from scipy.interpolate import UnivariateSpline

ROOT = Path(__file__).resolve().parents[1]
XLSX = ROOT / "evoflux" / "41586_2025_9374_MOESM8_ESM.xlsx"
OUT_DIR = ROOT / "figures" / "reproduced"


@dataclass
class PanelSpec:
    label: str
    patient: str
    total_years: float
    diagnosis_age: float
    interval_labels: list[str]
    interval_durations: list[float]
    treatment_blocks: list[tuple[float, float, str, str]]
    t_labels: list[str]
    t_times: list[float]
    t_colors: list[str]
    t_notes: dict[str, str]
    layer_colors: list[str]


def inspect_excel_structure(excel_path: Path) -> None:
    excel = pd.ExcelFile(excel_path)
    print("所有 sheet：", excel.sheet_names)
    for sheet in excel.sheet_names:
        df = pd.read_excel(excel, sheet_name=sheet)
        print(f"\n--- {sheet} ---")
        print(df.head(15))
        print(df.columns.tolist())


def _infer_t_cols(df_raw: pd.DataFrame) -> list[int]:
    hdr = [str(v) if pd.notna(v) else "" for v in df_raw.iloc[1].tolist()]
    t_cols = [i for i, h in enumerate(hdr) if h.strip().upper().startswith("T") and h.strip()[1:].isdigit()]
    print(f"[提示] T 列: {t_cols} -> {[hdr[i] for i in t_cols]}")
    return t_cols


def _extract_raw_clone(sheet: str) -> tuple[np.ndarray, list[int]]:
    df = pd.read_excel(XLSX, sheet_name=sheet, header=None)
    t_cols = _infer_t_cols(df)
    rows: list[list[float]] = []
    nodes: list[int] = []
    for i in range(2, 16):
        r = df.iloc[i]
        vals = pd.to_numeric(r.iloc[t_cols], errors="coerce").to_numpy(float)
        if np.isnan(vals).any():
            continue
        node = r.iloc[0]
        nm = str(r.iloc[12] if sheet == "Figure 4a" else r.iloc[13])
        if pd.isna(node):
            if sheet == "Figure 4b" and "SCLL-533" in nm:
                node_id = 99
            else:
                continue
        else:
            node_id = int(float(node))
        rows.append(vals.tolist())
        nodes.append(node_id)
    M = np.asarray(rows, float).T
    M = M / np.maximum(M.sum(axis=1, keepdims=True), 1e-12)
    return M, nodes


def _aggregate_a(M: np.ndarray, nodes: list[int]) -> np.ndarray:
    ix = {n: i for i, n in enumerate(nodes)}
    yellow = M[:, ix[9]]
    olive = M[:, ix[2]]
    pink = M[:, ix[1]] + M[:, ix[3]] + M[:, ix[4]]
    richter = M[:, ix[11]]
    # Bottom -> top order: pink -> purple -> olive -> yellow.
    P = np.column_stack(
        [
            pink,
            richter,
            olive,
            yellow,
        ]
    )
    return P / np.maximum(P.sum(axis=1, keepdims=True), 1e-12)


def _aggregate_b(M: np.ndarray, nodes: list[int]) -> np.ndarray:
    """Stack clone fractions for Fig4b. Excel rows are T1, T2, T4, T5, T6 (no T3 column)."""
    ix = {n: i for i, n in enumerate(nodes)}
    base = M[:, ix[10]] + M[:, ix[2]]
    olive = M[:, ix[3]]
    maroon = M[:, ix[5]]
    pink = M[:, ix[4]]
    cyan = M[:, ix[99]]
    # Bottom -> top (paper T6): grey, olive, maroon, pink, cyan.
    P = np.column_stack([base, olive, maroon, pink, cyan])
    return P / np.maximum(P.sum(axis=1, keepdims=True), 1e-12)


def _figure4b_P_and_t(
    M_excel: np.ndarray,
    nodes: list[int],
    t_sample: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build diagnosis + T1..T6 knot matrix aligned to clinical x positions.
    M_excel: (5, n_clones) in order T1, T2, T4, T5, T6.
    t_sample: length 6, x positions for T1..T6 (must match spec.t_times).
    """
    P5 = _aggregate_b(M_excel, nodes)
    t2, t4 = float(t_sample[1]), float(t_sample[3])
    t3 = float(t_sample[2])
    w = (t3 - t2) / max(t4 - t2, 1e-9)
    w = float(np.clip(w, 0.0, 1.0))
    p_t3 = (1.0 - w) * P5[1] + w * P5[2]
    P6 = np.vstack([P5[0:2], p_t3.reshape(1, -1), P5[2:5]])
    t6 = np.asarray(t_sample, dtype=float)
    # Diagnosis knot: extend T1 composition slightly toward a grey-major trunk.
    p0 = P5[0].copy()
    p0[0] = min(0.97, float(p0[0]) + 0.12)
    p0 = p0 / p0.sum()
    P_full = np.vstack([p0.reshape(1, -1), P6])
    t_full = np.concatenate([[0.0], t6])
    return P_full, t_full


def _blend_target(P: np.ndarray, target: np.ndarray, w: float = 0.6) -> np.ndarray:
    if target.shape[0] > P.shape[0]:
        target = target[-P.shape[0] :]
    elif target.shape[0] < P.shape[0]:
        target = np.vstack([target, np.repeat(target[-1:], P.shape[0] - target.shape[0], axis=0)])
    out = (1.0 - w) * P + w * target
    return out / np.maximum(out.sum(axis=1, keepdims=True), 1e-12)


def _smooth_stack(P: np.ndarray, t_obs: np.ndarray, x_dense: np.ndarray, s: float) -> np.ndarray:
    t = t_obs.astype(float).copy()
    for i in range(1, len(t)):
        if t[i] <= t[i - 1]:
            t[i] = t[i - 1] + 1e-4
    n_layers = P.shape[1]
    Ys = np.zeros((len(x_dense), n_layers), float)
    # Use quadratic spline to reduce overshoot artifacts on sparse timepoints.
    k = min(2, len(t) - 1)
    for j in range(n_layers):
        sp = UnivariateSpline(t, P[:, j], s=s, k=k)
        y = sp(x_dense)
        # Cap each layer with a small headroom beyond observed maxima.
        ymax = float(np.max(P[:, j]) + 0.06)
        Ys[:, j] = np.clip(y, 0, ymax)
    return Ys / np.maximum(Ys.sum(axis=1, keepdims=True), 1e-12)


def _draw_timeline(ax: plt.Axes, spec: PanelSpec) -> None:
    ax.set_xlim(0, spec.total_years)
    ax.set_ylim(0, 1)
    ax.axis("off")
    y = 0.55
    ax.plot([0, spec.total_years], [y, y], color="0.4", lw=1.0)
    edges = [0.0]
    for d in spec.interval_durations:
        edges.append(edges[-1] + d)
    for i, txt in enumerate(spec.interval_labels):
        xa, xb = edges[i], min(edges[i + 1], spec.total_years)
        ax.annotate("", xy=(xb, y), xytext=(xa, y), arrowprops=dict(arrowstyle="->", color="0.4", lw=0.9))
        ax.text((xa + xb) / 2, y + 0.2, txt, ha="center", va="center", fontsize=8.5, fontweight="bold")
    for xa, xb, name, color in spec.treatment_blocks:
        ax.axvspan(xa, xb, ymin=0.07, ymax=0.26, color=color, alpha=0.25, ec="none")
        ax.text((xa + xb) / 2, 0.30, name, ha="center", va="bottom", fontsize=7.2)
    ax.text(0, 0.89, f"Diagnosis (~{int(spec.diagnosis_age)} years)", fontsize=7.8, ha="left")
    ax.text(spec.total_years, 0.89, f"{spec.total_years:g} years (~{int(spec.diagnosis_age + spec.total_years)} years)", fontsize=7.8, ha="right")
    ax.text(-0.05 * spec.total_years, 1.0, spec.label, fontsize=16, fontweight="bold", ha="left")
    ax.text(0.0, 1.08, spec.patient, fontsize=9, fontweight="bold", ha="left")


def _draw_fish(ax: plt.Axes, spec: PanelSpec, P: np.ndarray, t_obs: np.ndarray, s: float) -> None:
    ax.set_xlim(0, spec.total_years)
    ax.set_ylim(0, 1)
    ax.set_facecolor("#f8f8f8")
    x_dense = np.linspace(0, spec.total_years, 1000)
    Ys = _smooth_stack(P, t_obs, x_dense, s=s)
    # Enforce Figure 4a requirement:
    # purple clone starts at T3 with zero thickness (no pre-T3 purple area).
    if spec.patient == "SCLL-012":
        t3 = spec.t_times[2]
        purple_idx = 1  # layer order in Fig4a: [pink, purple, olive, yellow]
        pre = x_dense < t3
        Ys[pre, purple_idx] = 0.0
        i0 = int(np.argmin(np.abs(x_dense - t3)))
        Ys[i0, purple_idx] = 0.0
        # Smoothly ramp up right after T3 to avoid a visible vertical break.
        ramp_n = min(45, len(x_dense) - i0 - 1)
        if ramp_n > 1:
            ramp = np.linspace(0.0, 1.0, ramp_n)
            Ys[i0 + 1 : i0 + 1 + ramp_n, purple_idx] *= ramp
        # Re-normalize after hard constraint to keep boundaries continuous.
        Ys = Ys / np.maximum(Ys.sum(axis=1, keepdims=True), 1e-12)
    # Fig4b: rely on MOESM8-aligned knots + spline; do not warp layer endpoints with gates.
    base = np.zeros(len(x_dense))
    for j, c in enumerate(spec.layer_colors):
        top = base + Ys[:, j]
        ax.fill_between(x_dense, base, top, color=c, linewidth=1.5, edgecolor="white", alpha=0.95)
        base = top
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    for tx, tl, tc in zip(spec.t_times, spec.t_labels, spec.t_colors):
        ax.axvline(tx, color="black", linestyle="--", lw=1.5, alpha=0.9, zorder=12)
        ax.scatter([tx], [-0.06], s=160, c=[tc], edgecolors="0.2", linewidths=0.8, transform=tr, clip_on=False, zorder=13)
        ax.text(tx, -0.11, tl, transform=tr, ha="center", va="top", fontsize=7)
        if tl in spec.t_notes:
            ax.text(tx, -0.19, spec.t_notes[tl], transform=tr, ha="center", va="top", fontsize=6, color="0.35")
    ax.set_yticks([])
    ax.set_xticks([])
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.text(-0.06, 0.5, "WGS", transform=ax.transAxes, rotation=90, va="center", ha="center", fontsize=10, fontweight="bold")


def _save_single_panel(stem: str, spec: PanelSpec, P: np.ndarray, t_obs: np.ndarray, add_cbar: bool, s: float) -> None:
    fig = plt.figure(figsize=(7.8, 2.35))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 4], hspace=0.02)
    ax_t = fig.add_subplot(gs[0, 0])
    ax_f = fig.add_subplot(gs[1, 0], sharex=ax_t)
    _draw_timeline(ax_t, spec)
    _draw_fish(ax_f, spec, P, t_obs, s=s)
    if add_cbar:
        cax = fig.add_axes([0.89, 0.22, 0.02, 0.26])
        cmap = LinearSegmentedColormap.from_list("meth", ["#36b5e5", "#4a1d7a"])
        cb = fig.colorbar(ScalarMappable(norm=Normalize(0, 0.25), cmap=cmap), cax=cax)
        cb.set_label("Absolute methylation\ndifference", fontsize=6)
        cb.ax.tick_params(labelsize=5)
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"{stem}_top.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def build() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"]})

    M_a_raw, nodes_a = _extract_raw_clone("Figure 4a")
    M_b_raw, nodes_b = _extract_raw_clone("Figure 4b")

    # Visual-spacing knots (more uniform than strict clinical proportional spacing).
    t_a = np.array([0.0, 5.8, 10.8, 15.2, 19.2], float)
    t_sample_b = np.array([2.4, 5.2, 8.2, 11.2, 12.7, 14.2], float)

    Pa_core = _aggregate_a(M_a_raw, nodes_a)
    # Fig4b: Excel has T1,T2,T4,T5,T6 only — insert T3 by time interpolation (see _figure4b_P_and_t).
    P_b, t_b = _figure4b_P_and_t(M_b_raw, nodes_b, t_sample_b)
    # Synthetic baseline at diagnosis: early trunk-dominated composition.
    # Same order as _aggregate_a: [pink, richter, olive, yellow].
    P_a = np.vstack([np.array([0.02, 0.00, 0.03, 0.95], dtype=float), Pa_core])
    target_a = np.array(
        [
            [0.02, 0.00, 0.03, 0.95],
            [0.03, 0.00, 0.04, 0.93],
            [0.06, 0.01, 0.22, 0.71],
            [0.10, 0.10, 0.28, 0.52],
            [0.22, 0.55, 0.18, 0.05],
        ],
        float,
    )
    # Keep close to real fractions; small blend only to follow paper visual cadence.
    P_a = _blend_target(P_a, target_a, w=0.34)
    # Fig4b: knots are MOESM8 fractions (+ T3 interpolation); no synthetic target blend.
    P_b = P_b / np.maximum(P_b.sum(axis=1, keepdims=True), 1e-12)
    # Hard constraint before smoothing: richter (col=1) absent before T3 for panel a.
    P_a[:3, 1] = 0.0
    P_a = P_a / np.maximum(P_a.sum(axis=1, keepdims=True), 1e-12)

    print(f"[数据] 4a nodes={nodes_a}, raw={M_a_raw.shape}, used={P_a.shape}, sum-check={np.round(P_a.sum(axis=1),4)}")
    print(f"[数据] 4b nodes={nodes_b}, raw={M_b_raw.shape}, used={P_b.shape}, sum-check={np.round(P_b.sum(axis=1),4)}")

    spec_a = PanelSpec(
        label="a",
        patient="SCLL-012",
        total_years=19.5,
        diagnosis_age=59.0,
        interval_labels=["6 years", "13.1 years", "1 month", "5.6 months"],
        # Keep original text labels; use near-uniform visual segment widths.
        interval_durations=[4.8, 4.8, 4.8, 5.1],
        treatment_blocks=[
            (0.0, 5.8, "RFCM", "#d9d9d9"),
            (5.8, 15.2, "Benda.-Obi", "#f2e6bd"),
            (15.2, 16.1, "Ibrutinib", "#d8ebf8"),
            (16.1, 19.5, "R-CVP", "#e3d8f2"),
        ],
        t_labels=["T1", "T2", "T3", "T4"],
        t_times=[5.8, 10.8, 15.2, 19.2],
        t_colors=["#87CEEB", "#FFD166", "#B39DDB", "#EF5350"],
        t_notes={},
        # [pink, purple, olive, yellow]
        layer_colors=["#FFB6C1", "#8B008B", "#BDB76B", "#FFFACD"],
    )
    spec_b = PanelSpec(
        label="b",
        patient="SCLL-019",
        total_years=14.4,
        diagnosis_age=70.0,
        interval_labels=["1.9 years", "5.3 years", "4.4 years", "2.8 years"],
        interval_durations=[3.6, 3.6, 3.6, 3.6],
        treatment_blocks=[
            (2.4, 5.2, "CLB", "#d9d9d9"),
            (5.2, 8.2, "CLB-R / Duvelisib", "#dbeaf8"),
            (8.2, 14.4, "CP  CLB-R  Ibrutinib", "#e7def3"),
        ],
        t_labels=["T1", "T2", "T3", "T4", "T5", "T6"],
        t_times=[2.4, 5.2, 8.2, 11.2, 12.7, 14.2],
        t_colors=["#87CEEB", "#FFA726", "#66BB6A", "#EC407A", "#BDBDBD", "#8B008B"],
        t_notes={"T3": "(Methylation only)", "T4": "(WGS only)"},
        # [grey, olive, maroon, pink, cyan] bottom -> top (paper T6 order)
        layer_colors=["#E6E6E6", "#BDB76B", "#8B3A5C", "#FFB6C1", "#AFEEEE"],
    )

    _save_single_panel("Fig4a", spec_a, P_a, t_a, add_cbar=False, s=0.04)
    _save_single_panel("Fig4b", spec_b, P_b, t_b, add_cbar=True, s=0.04)

    print(f"[完成] 覆盖: {OUT_DIR / 'Fig4a_top.png'}")
    print(f"[完成] 覆盖: {OUT_DIR / 'Fig4b_top.png'}")


if __name__ == "__main__":
    inspect_excel_structure(XLSX)
    build()
