from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
MOESM3 = ROOT / "evoflux" / "41586_2025_9374_MOESM3_ESM.xlsx"
OUT_DIR = ROOT / "figures" / "reproduced" / "moesm3_methods_processed"
REFERENCE_BETA = ROOT / "data" / "beta_fcpgs.csv"


TABLES = {
    "supplementary_table_2": "step1_2_sample_metadata",
    "supplementary_table_3": "step1_2_qc_deconv",
    "supplementary_table_4": "step3_fcpg_annotation",
    "supplementary_table_5": "step3_fcpg_beta_2204",
    "supplementary_table_6": "step4_cna_cll",
    "supplementary_table_7": "step4_cna_mcl",
    "supplementary_table_12": "step10_11_evoflux_inference",
    "supplementary_table_18": "step12_longitudinal_rt_samples",
}


def _read_sheet(sheet: str) -> pd.DataFrame:
    # MOESM3 uses 2 descriptive rows, then data header.
    df = pd.read_excel(MOESM3, sheet_name=sheet, header=2)
    df = df.dropna(axis=0, how="all").dropna(axis=1, how="all")
    return df


def _write_csv(df: pd.DataFrame, name: str) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / f"{name}.csv"
    df.to_csv(out, index=False)
    return out


def _process_fcpg_beta(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    beta = df.rename(columns={"fCpGs": "cpg_id"}).copy()
    beta = beta[beta["cpg_id"].notna()].copy()
    beta["cpg_id"] = beta["cpg_id"].astype(str)
    beta = beta.set_index("cpg_id")
    beta = beta.apply(pd.to_numeric, errors="coerce")

    summary: dict[str, object] = {
        "rows": int(beta.shape[0]),
        "cols": int(beta.shape[1]),
    }

    if REFERENCE_BETA.exists():
        ref = pd.read_csv(REFERENCE_BETA, index_col=0)
        common_rows = beta.index.intersection(ref.index)
        common_cols = beta.columns.intersection(ref.columns)
        summary["reference_rows"] = int(ref.shape[0])
        summary["reference_cols"] = int(ref.shape[1])
        summary["common_rows"] = int(len(common_rows))
        summary["common_cols"] = int(len(common_cols))
        if len(common_rows) and len(common_cols):
            diff = (beta.loc[common_rows, common_cols] - ref.loc[common_rows, common_cols]).abs()
            summary["max_abs_diff_vs_reference"] = float(diff.max().max())
            summary["mean_abs_diff_vs_reference"] = float(diff.stack().mean())
    return beta, summary


def main() -> None:
    if not MOESM3.exists():
        raise FileNotFoundError(f"MOESM3 not found: {MOESM3}")

    outputs: dict[str, str] = {}
    checks: dict[str, object] = {}
    loaded: dict[str, pd.DataFrame] = {}

    for sheet, tag in TABLES.items():
        df = _read_sheet(sheet)
        loaded[sheet] = df
        out = _write_csv(df, tag)
        outputs[sheet] = str(out)
        checks[f"{sheet}_shape"] = [int(df.shape[0]), int(df.shape[1])]

    # Structured derived outputs used by downstream reproduction steps.
    st2 = loaded["supplementary_table_2"].rename(columns={"participant_id_anonymous": "PARTICIPANT_ID_ANONYMOUS"})
    st3 = loaded["supplementary_table_3"].rename(columns={"participant_id_anonymous": "PARTICIPANT_ID_ANONYMOUS"})
    st12 = loaded["supplementary_table_12"]

    merge_cols_2 = [
        "PARTICIPANT_ID_ANONYMOUS",
        "sample_id",
        "cell_type",
        "cell_type_annot_1",
        "cell_type_annot_2",
        "cell_type_annot_3",
        "disease_subtype",
        "disease_subtype_epigenetics",
    ]
    merge_cols_3 = [
        "PARTICIPANT_ID_ANONYMOUS",
        "sample_group_analysis_normalization",
        "ssnob_normalization_batch",
        "purity_tumor_consensus",
        "age_sampling",
    ]

    st12_merged = (
        st12.merge(st2[[c for c in merge_cols_2 if c in st2.columns]], on="PARTICIPANT_ID_ANONYMOUS", how="left")
        .merge(st3[[c for c in merge_cols_3 if c in st3.columns]], on="PARTICIPANT_ID_ANONYMOUS", how="left")
    )
    outputs["derived_step10_11_inference_with_metadata"] = str(
        _write_csv(st12_merged, "derived_step10_11_inference_with_metadata")
    )
    checks["derived_step10_11_inference_with_metadata_shape"] = [int(st12_merged.shape[0]), int(st12_merged.shape[1])]

    beta, beta_summary = _process_fcpg_beta(loaded["supplementary_table_5"])
    beta_out = OUT_DIR / "derived_step3_fcpg_beta_matrix.csv"
    beta.to_csv(beta_out)
    outputs["derived_step3_fcpg_beta_matrix"] = str(beta_out)
    checks["derived_step3_fcpg_beta_matrix"] = beta_summary

    # Helpful extraction for Fig4 / PISCA longitudinal sample IDs.
    selected = ["SCLL-545", "SCLL-546", "SCLL-547", "SCLL-548", "SCLL-493", "SCLL-494", "SCLL-531", "SCLL-532", "SCLL-533"]
    present = [c for c in selected if c in beta.columns]
    beta_pisca = beta[present].copy()
    pisca_out = OUT_DIR / "derived_step12_pisca_samples_beta.csv"
    beta_pisca.to_csv(pisca_out)
    outputs["derived_step12_pisca_samples_beta"] = str(pisca_out)
    checks["derived_step12_pisca_samples_beta_shape"] = [int(beta_pisca.shape[0]), int(beta_pisca.shape[1])]

    report = {
        "moesm3_path": str(MOESM3),
        "outputs": outputs,
        "checks": checks,
        "notes": [
            "Tables are extracted from MOESM3 using header row index 2.",
            "The fCpG beta matrix is compared against data/beta_fcpgs.csv when available.",
            "Derived PISCA sample matrix is included for longitudinal phylogenetic workflows.",
        ],
    }
    report_path = OUT_DIR / "processing_report.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"WROTE {report_path}")


if __name__ == "__main__":
    main()

