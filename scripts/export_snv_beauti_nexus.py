"""
Export MOESM3-derived binary SNV matrix + tip ages to NEXUS (DATATYPE=BINARY) for BEAST 1.10 Beauti.

The paper uses BEAST 1.10 for SNV phylogeny; PISCA is not used here (see PISCA README).

This does not run BEAST. Provide ages that match your analysis (e.g. from clinical metadata).

Usage:
  python scripts/export_snv_beauti_nexus.py \\
    --matrix figures/reproduced/snv_beast10_prep/case_1_binary_matrix.csv \\
    --ages-csv path/to/ages.csv \\
    --out figures/reproduced/snv_beast10_prep/case_1_for_beauti.nexus

ages.csv columns: Name (sample ID matching matrix rows), AGE_SAMPLING (years).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser(description="Binary SNV matrix → NEXUS for Beauti / BEAST 1.10.")
    ap.add_argument("--matrix", type=Path, required=True, help="Rows = samples, cols = sites, values 0/1")
    ap.add_argument("--ages-csv", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    mat = pd.read_csv(args.matrix, index_col=0)
    ages = pd.read_csv(args.ages_csv)
    if "Name" not in ages.columns or "AGE_SAMPLING" not in ages.columns:
        raise SystemExit("ages.csv must have columns: Name, AGE_SAMPLING")

    samples = [str(x) for x in mat.index.tolist()]
    age_map = ages.set_index("Name")["AGE_SAMPLING"].to_dict()
    missing = [s for s in samples if s not in age_map]
    if missing:
        raise SystemExit(f"ages.csv missing samples: {missing}")

    ntax = len(samples)
    nchar = mat.shape[1]
    lines = [
        "#NEXUS",
        "BEGIN DATA;",
        f"\tDIMENSIONS NTAX={ntax} NCHAR={nchar};",
        "\tFORMAT DATATYPE=BINARY MISSING=? GAP=-;",
        "\tMATRIX",
    ]
    for sid in samples:
        row = mat.loc[sid].astype(int)
        seq = "".join(str(int(x)) for x in row.values)
        taxon_id = sid.replace(" ", "_").replace("-", "_")
        lines.append(f"\t\t{taxon_id}\t{seq}")
    lines.append("\t;")
    lines.append("END;")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines), encoding="utf-8")
    sidecar = args.out.with_suffix(".tip_ages.tsv")
    pd.DataFrame(
        [{"Name": s, "AGE_SAMPLING": float(age_map[s])} for s in samples]
    ).to_csv(sidecar, sep="\t", index=False)
    print(f"Wrote {args.out} (NTAX={ntax}, NCHAR={nchar})")
    print(f"Wrote tip ages: {sidecar}")
    print("Import the NEXUS in BEAUTi 1.10.x, set tip dates from the .tip_ages.tsv, then choose clock + tree prior.")


if __name__ == "__main__":
    main()
