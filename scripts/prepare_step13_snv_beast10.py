from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
_SCRIPTS = ROOT / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

MOESM3 = ROOT / "evoflux" / "41586_2025_9374_MOESM3_ESM.xlsx"
OUT = ROOT / "figures" / "reproduced" / "snv_beast10_prep"


def _check_beast10() -> dict[str, object]:
    candidates = [
        ROOT / "tools" / "BEAST v1.10.0",
        ROOT / "tools" / "BEAST v1.10.4",
        ROOT / "tools" / "BEAST v1.10.5",
        Path(r"D:\tools\BEASTv1.10.5"),
        Path(r"C:\tools\BEASTv1.10.5"),
    ]
    found = None
    for c in candidates:
        if (c / "bin" / "beast.cmd").exists() and (c / "lib" / "beast.jar").exists():
            found = c
            break
    return {
        "beast10_found": str(found) if found else None,
        "checked_candidates": [str(x) for x in candidates],
    }


def _site_id(df: pd.DataFrame) -> pd.Series:
    return (
        df["CHROM"].astype(str)
        + ":"
        + df["POSITION"].astype(int).astype(str)
        + "_"
        + df["REF"].astype(str)
        + ">"
        + df["ALT"].astype(str)
    )


def _build_case_matrix(case_df: pd.DataFrame) -> pd.DataFrame:
    samples = sorted(case_df["SAMPLE"].dropna().astype(str).unique().tolist())
    case_df = case_df.copy()
    case_df["site_id"] = _site_id(case_df)
    mut = case_df[["SAMPLE", "site_id"]].drop_duplicates()
    mut["present"] = 1
    mat = mut.pivot(index="SAMPLE", columns="site_id", values="present").fillna(0).astype(int)
    # Ensure stable row order
    mat = mat.loc[samples]
    # Remove invariant sites (all 1 here cannot happen; all 0 impossible due pivot),
    # but keep only polymorphic sites with at least one 1 and one 0 across samples.
    keep = (mat.sum(axis=0) > 0) & (mat.sum(axis=0) < mat.shape[0])
    return mat.loc[:, keep]


def _sbs1_feasibility(columns: list[str]) -> dict[str, object]:
    """COSMIC SBS1 needs trinucleotide context; MOESM3 table_15 often lacks it."""
    cset = {str(x).upper() for x in columns}
    hints = ("TRINUC", "TRIPLE", "CONTEXT", "5'", "3'", "STRAND")
    has_ctx = any(any(h in str(col).upper() for h in hints) for col in columns)
    return {
        "supplementary_table_15_has_trinucleotide_columns": bool(has_ctx),
        "sbs1_filtering_feasible_from_moesm3_alone": bool(has_ctx),
        "note": (
            "SBS1 (clock-like) filtering typically uses reference-genome trinucleotide context "
            "and/or signature assignment (e.g. SigProfiler). MOESM3 supplementary_table_15 "
            "as shipped here only has REF/ALT/position — insufficient for strict SBS1 without "
            "lifting to hg38 or author-provided signature calls."
        ),
    }


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Extract WGS SNVs from MOESM3 and build 0/1 site matrices for BEAST 1.10 prep."
    )
    ap.add_argument(
        "--hg38-fasta",
        type=Path,
        default=None,
        help="Optional hg38.fa (indexed). Adds ref trinucleotide columns via scripts/annotate_snv_trinucleotide.py",
    )
    ap.add_argument(
        "--chrom-prefix",
        default="chr",
        help="Chromosome prefix for FASTA (default chr when Excel uses numeric CHROM).",
    )
    ap.add_argument(
        "--filter-cpg-c-to-t",
        action="store_true",
        help=(
            "Restrict to C>T at CpG in hg38 (requires --hg38-fasta). "
            "This is a clock-like subset proxy — not full COSMIC SBS1 / SigProfiler assignment."
        ),
    )
    args = ap.parse_args()

    if args.filter_cpg_c_to_t and not args.hg38_fasta:
        raise SystemExit("--filter-cpg-c-to-t requires --hg38-fasta")

    OUT.mkdir(parents=True, exist_ok=True)

    if not MOESM3.exists():
        raise FileNotFoundError(f"MOESM3 not found: {MOESM3}")

    # WGS SNVs + WGS mode annotations
    t15 = pd.read_excel(MOESM3, sheet_name="supplementary_table_15", header=2)
    t16 = pd.read_excel(MOESM3, sheet_name="supplementary_table_16", header=2)
    sbs1_info = _sbs1_feasibility(list(t15.columns))

    # Keep SNVs only; this is still not SBS1-filtered.
    wgs = t15[t15["TYPE"].astype(str).str.upper() == "SNV"].copy()
    wgs = wgs.dropna(subset=["CASE", "SAMPLE", "CHROM", "POSITION", "REF", "ALT"])

    annotation_path: str | None = None
    sbs1_filter_note = None
    if args.hg38_fasta:
        from annotate_snv_trinucleotide import annotate_dataframe

        wgs = annotate_dataframe(wgs, args.hg38_fasta, chrom_prefix=args.chrom_prefix)
        annotation_path = str(OUT / "supplementary_table_15_wgs_snv_annotated.csv")
        wgs.to_csv(annotation_path, index=False)
        sbs1_info = {
            **sbs1_info,
            "hg38_annotation_applied": True,
            "annotated_csv": annotation_path,
            "filter_cpg_c_to_t_applied": bool(args.filter_cpg_c_to_t),
        }
        if args.filter_cpg_c_to_t:
            n0 = len(wgs)
            wgs = wgs[wgs["is_cpg_c_to_t"]].copy()
            sbs1_filter_note = f"Filtered SNVs to is_cpg_c_to_t: {n0} -> {len(wgs)} rows"
            sbs1_info["cpg_c_to_t_filter_note"] = sbs1_filter_note

    # Candidate cases for phylogeny: at least 2 longitudinal samples.
    case_n = wgs.groupby("CASE")["SAMPLE"].nunique().sort_values(ascending=False)
    multi_cases = case_n[case_n >= 2]

    outputs: dict[str, object] = {
        "wgs_snv_rows": int(wgs.shape[0]),
        "n_unique_cases": int(wgs["CASE"].nunique()),
        "n_cases_with_ge2_samples": int(multi_cases.shape[0]),
        "top_cases_by_sample_count": {str(k): int(v) for k, v in multi_cases.head(10).items()},
    }

    # Build 0/1 matrices for top 5 multi-sample cases.
    case_files: dict[str, str] = {}
    for case_id in multi_cases.head(5).index.tolist():
        case_df = wgs[wgs["CASE"] == case_id]
        mat = _build_case_matrix(case_df)
        out_csv = OUT / f"case_{case_id}_binary_matrix.csv"
        mat.to_csv(out_csv)
        case_files[str(case_id)] = str(out_csv)
        outputs[f"case_{case_id}_shape"] = [int(mat.shape[0]), int(mat.shape[1])]

    # Persist raw extracts useful for downstream filtering / BEAST.
    wgs_out = OUT / "supplementary_table_15_wgs_snv.csv"
    t16_out = OUT / "supplementary_table_16_wgs_modes.csv"
    wgs.to_csv(wgs_out, index=False)
    t16.to_csv(t16_out, index=False)

    beast10 = _check_beast10()
    sbs1_done = bool(args.hg38_fasta and args.filter_cpg_c_to_t)
    next_req = [
        "Optional: scripts/export_snv_beauti_nexus.py --matrix case_*_binary_matrix.csv --ages-csv ... --out ..._for_beauti.nexus (NEXUS + .tip_ages.tsv for BEAUTi 1.10 import).",
        "Install BEAST 1.10.x and build XML (e.g. Binary + strict clock + constant coalescent + tip dates) in Beauti or a validated template.",
        "Run MCMC and assess ESS in Tracer.",
    ]
    if not args.hg38_fasta:
        next_req.insert(
            0,
            "Optional: pass --hg38-fasta to add trinucleotide context; use --filter-cpg-c-to-t for a CpG C>T clock-like subset (not full SigProfiler SBS1).",
        )

    report = {
        "moesm3": str(MOESM3),
        "outputs": {
            "wgs_snv_table": str(wgs_out),
            "wgs_mode_table": str(t16_out),
            "case_binary_matrices": case_files,
            "wgs_snv_annotated": annotation_path,
        },
        "summary": outputs,
        "beast10_preflight": beast10,
        "status": {
            "sbs1_proxy_filtering_done": sbs1_done,
            "beast10_xml_built": False,
            "beast10_run_done": False,
        },
        "sbs1_feasibility": sbs1_info,
        "next_required": next_req,
    }
    rep = OUT / "step13_preflight_report.json"
    rep.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"WROTE {rep}")


if __name__ == "__main__":
    main()

