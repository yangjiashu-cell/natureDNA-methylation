"""
Add reference trinucleotide context to MOESM3-style WGS SNV rows (hg38).

Requires: pyfaidx (`pip install pyfaidx`) and an indexed hg38 FASTA (e.g. hg38.fa + hg38.fa.fai).

This does NOT run SigProfiler / full COSMIC SBS1 assignment. It only adds:
  - ref_trinuc: 3bp on forward strand at (POS-1, POS, POS+1) with REF at the middle
  - is_c_to_t: REF==C and ALT==T
  - is_cpg_c_to_t: C>T at a CpG dinucleotide (middle C and 3' G in reference)

Use the latter as a *biologically motivated clock-like subset* when the supplementary table
lacks trinucleotide columns — not a substitute for signature deconvolution in the paper.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def _norm_chrom(c: str | int, prefix: str) -> str:
    s = str(c).strip()
    if s.lower().startswith("chr"):
        return s
    return f"{prefix}{s}"


def annotate_dataframe(df: pd.DataFrame, fasta_path: Path, *, chrom_prefix: str = "chr") -> pd.DataFrame:
    try:
        from pyfaidx import Fasta
    except ImportError as e:
        raise RuntimeError("Install pyfaidx to use --hg38-fasta: pip install pyfaidx") from e

    fa = Fasta(str(fasta_path))
    trinucs: list[str] = []
    c_to_t: list[bool] = []
    cpg_ct: list[bool] = []

    for _, row in df.iterrows():
        chrom = _norm_chrom(row["CHROM"], chrom_prefix)
        pos = int(row["POSITION"])
        ref = str(row["REF"]).upper()
        alt = str(row["ALT"]).upper()
        if chrom not in fa:
            trinucs.append("")
            c_to_t.append(False)
            cpg_ct.append(False)
            continue
        # pyfaidx: 1-based inclusive slice [pos-1 : pos+1] => three bases around SNV
        try:
            # pyfaidx: 1-based inclusive slice → bases at POS-1, POS, POS+1
            tri = str(fa[chrom][pos - 1 : pos + 1]).upper()
        except (KeyError, ValueError):
            trinucs.append("")
            c_to_t.append(False)
            cpg_ct.append(False)
            continue
        if len(tri) != 3:
            trinucs.append("")
            c_to_t.append(False)
            cpg_ct.append(False)
            continue
        tri_ok = tri[1] == ref
        trinucs.append(tri if tri_ok else tri + "?")
        ic = ref == "C" and alt == "T"
        c_to_t.append(ic)
        # CpG on forward: mutated C with G immediately downstream
        cpg_ct.append(ic and tri[1] == "C" and tri[2] == "G")

    out = df.copy()
    out["ref_trinuc"] = trinucs
    out["is_c_to_t"] = c_to_t
    out["is_cpg_c_to_t"] = cpg_ct
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Annotate WGS SNV table with hg38 trinucleotide context.")
    ap.add_argument("--wgs-csv", type=Path, required=True, help="e.g. supplementary_table_15_wgs_snv.csv")
    ap.add_argument("--hg38-fasta", type=Path, required=True, help="Path to hg38.fa (must be indexed).")
    ap.add_argument("--chrom-prefix", default="chr", help="Prefix if Excel uses numeric CHROM (default chr).")
    ap.add_argument("--out-csv", type=Path, required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.wgs_csv)
    ann = annotate_dataframe(df, args.hg38_fasta, chrom_prefix=args.chrom_prefix)
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    ann.to_csv(args.out_csv, index=False)
    print(f"Wrote {args.out_csv} ({ann.shape[0]} rows)")


if __name__ == "__main__":
    main()
