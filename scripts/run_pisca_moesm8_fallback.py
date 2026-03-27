"""
MOESM8 Fig.4a/4b fallback when Zenodo harmonised matrix is unavailable.

Discretisation uses fixed β thresholds (0.2 / 0.8), not the paper’s Stan Beta mixture.
For Methods-aligned methylation phylogeny, prefer scripts/reproduce_pisca_zenodo_fig4.py
and scripts/run_fig4_pisca_paper_methods_pipeline.py with data/beta_fcpgs.csv.
"""
from __future__ import annotations

import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MOESM8 = ROOT / "evoflux" / "41586_2025_9374_MOESM8_ESM.xlsx"
BEAST_CMD = ROOT / "tools" / "BEAST v1.8.4" / "bin" / "beast.cmd"
TEMPLATE_XML = ROOT / "tools" / "PISCAv1.1" / "examples" / "biallelicBinary4Params.xml"
OUT_DIR = ROOT / "figures" / "reproduced" / "pisca_moesm8_fallback"


def _parse_sample_pair(text: str) -> tuple[str, str]:
    m1 = re.search(r"sample\.t1:([A-Za-z0-9_-]+)", text)
    m2 = re.search(r"sample\.t2:([A-Za-z0-9_-]+)", text)
    if not m1 or not m2:
        raise ValueError(f"Cannot parse sample pair: {text}")
    return m1.group(1), m2.group(1)


def _build_matrix(sheet_name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_excel(MOESM8, sheet_name=sheet_name, header=1)

    # Sample metadata block (age in years at sampling).
    meta = df[["Name", "AGE_SAMPLING"]].dropna(subset=["Name", "AGE_SAMPLING"]).copy()
    meta = meta[meta["Name"].astype(str).str.startswith("SCLL-")]
    meta = meta[~meta["Name"].astype(str).str.contains("newick-tree", case=False, na=False)]
    meta = meta.drop_duplicates(subset=["Name"]).reset_index(drop=True)

    # Pairwise methylation block.
    pairs = df[["CpG_name", "samples", "meth.T1", "meth.T2"]].dropna(subset=["CpG_name", "samples"]).copy()
    pairs["sample1"], pairs["sample2"] = zip(*pairs["samples"].map(_parse_sample_pair))

    long_rows: list[dict[str, object]] = []
    for _, r in pairs.iterrows():
        long_rows.append({"sample": r["sample1"], "CpG_name": r["CpG_name"], "beta": float(r["meth.T1"])})
        long_rows.append({"sample": r["sample2"], "CpG_name": r["CpG_name"], "beta": float(r["meth.T2"])})
    long_df = pd.DataFrame(long_rows)

    mat = long_df.pivot_table(index="sample", columns="CpG_name", values="beta", aggfunc="mean")
    mat = mat.dropna(axis=1, how="any")
    mat = mat.loc[[s for s in meta["Name"].tolist() if s in mat.index]]

    return mat, meta


def _beta_to_state(v: float) -> str:
    if v < 0.2:
        return "0"
    if v > 0.8:
        return "2"
    return "1"


def _replace_taxa_and_alignment(root: ET.Element, mat: pd.DataFrame, meta: pd.DataFrame) -> None:
    taxa = root.find("taxa")
    alignment = root.find("alignment")
    if taxa is None or alignment is None:
        raise RuntimeError("Template XML is missing <taxa> or <alignment> block")

    for child in list(taxa):
        taxa.remove(child)
    for child in list(alignment):
        alignment.remove(child)

    data_type_ref = ET.SubElement(alignment, "dataType")
    data_type_ref.set("idref", "biallelicBinary")

    for sample in mat.index:
        age = float(meta.loc[meta["Name"] == sample, "AGE_SAMPLING"].iloc[0])
        taxon = ET.SubElement(taxa, "taxon")
        taxon.set("id", sample.replace("-", "_"))
        date = ET.SubElement(taxon, "date")
        date.set("value", f"{age:.6f}")
        date.set("direction", "forwards")
        date.set("units", "years")

        seq = "".join(_beta_to_state(float(v)) for v in mat.loc[sample].values)
        sequence = ET.SubElement(alignment, "sequence")
        taxon_ref = ET.SubElement(sequence, "taxon")
        taxon_ref.set("idref", sample.replace("-", "_"))
        taxon_ref.tail = f"\n\t\t{seq}\n\t"


def _tune_mcmc_for_smoke_test(root: ET.Element, stem: str, n_taxa: int, max_age: float) -> None:
    mcmc = root.find("mcmc")
    if mcmc is None:
        raise RuntimeError("Template XML is missing <mcmc>")
    mcmc.set("chainLength", "50000")

    for p in root.findall(".//parameter"):
        pid = p.get("id")
        if pid == "luca_height":
            p.set("value", f"{max_age + 1.0:.6f}")
        elif pid == "luca_branch":
            p.set("value", "5")
            p.set("upper", "200")
            p.set("lower", "0.0")

    for u in root.findall(".//uniformPrior"):
        if u.find(".//parameter[@idref='luca_branch']") is not None:
            u.set("upper", "200")
            u.set("lower", "0.0")

    # Exchange operators are invalid for 2-taxon trees.
    if n_taxa < 3:
        operators = root.find("operators")
        if operators is not None:
            for child in list(operators):
                if child.tag in {
                    "wideExchange",
                    "narrowExchange",
                    "wilsonBalding",
                    "subtreeSlide",
                    "uniformOperator",
                    "upDownOperator",
                }:
                    operators.remove(child)

    for node in root.findall(".//log"):
        if node.get("id") == "screenLog":
            node.set("logEvery", "2000")
        if node.get("id") == "fileLog":
            node.set("logEvery", "2000")
            node.set("fileName", f"{stem}.log")
            node.set("overwrite", "true")
    for node in root.findall(".//logTree"):
        node.set("logEvery", "2000")
        node.set("fileName", f"{stem}.trees")


def build_case_xml(sheet_name: str, stem: str) -> Path:
    mat, meta = _build_matrix(sheet_name)
    if mat.empty or mat.shape[0] < 2 or mat.shape[1] < 100:
        raise RuntimeError(f"{sheet_name}: insufficient matrix after cleanup: {mat.shape}")

    tree = ET.parse(TEMPLATE_XML)
    root = tree.getroot()
    _replace_taxa_and_alignment(root, mat, meta)
    _tune_mcmc_for_smoke_test(root, stem, n_taxa=mat.shape[0], max_age=float(meta["AGE_SAMPLING"].max()))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    xml_path = OUT_DIR / f"{stem}.xml"
    tree.write(xml_path, encoding="utf-8", xml_declaration=True)

    # Also export matrix used in this fallback run for traceability.
    mat.to_csv(OUT_DIR / f"{stem}_matrix_beta.csv")
    return xml_path


def run_beast(xml_path: Path) -> int:
    cmd = [str(BEAST_CMD), "-overwrite", "-beagle_off", str(xml_path)]
    completed = subprocess.run(cmd, cwd=str(xml_path.parent), check=False)
    return completed.returncode


def main() -> None:
    jobs = [("Figure 4a", "fig4a_case12_moesm8_fallback"), ("Figure 4b", "fig4b_case19_moesm8_fallback")]
    for sheet, stem in jobs:
        xml_path = build_case_xml(sheet, stem)
        code = run_beast(xml_path)
        print(f"{sheet} -> {xml_path} | BEAST exit_code={code}")


if __name__ == "__main__":
    main()
