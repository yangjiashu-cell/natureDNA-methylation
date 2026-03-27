from __future__ import annotations

import argparse
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
NEXUS_IN = ROOT / "figures" / "reproduced" / "snv_beast_from_moesm3" / "alignment_binary.nexus"
AGES_CASE12 = ROOT / "figures" / "reproduced" / "pisca_zenodo_fig4" / "fig4a_case12_zenodo_ages.csv"
AGES_CASE19 = ROOT / "figures" / "reproduced" / "pisca_zenodo_fig4" / "fig4b_case19_zenodo_ages.csv"
OUT_DIR = ROOT / "figures" / "reproduced" / "snv_beast10_prep"


def _find_beast10_cmd() -> Path:
    p = ROOT / "tools" / "BEAST v1.10.4" / "bin" / "beast.cmd"
    if p.exists():
        return p
    raise FileNotFoundError(p)


def _parse_nexus(path: Path) -> tuple[list[str], dict[str, str]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    m = re.search(r"(?is)\bmatrix\b(.*?)\s*;", text)
    if not m:
        raise ValueError("No MATRIX block")
    taxa: list[str] = []
    seqs: dict[str, str] = {}
    for ln in m.group(1).strip().splitlines():
        parts = re.split(r"\s+", ln.strip())
        if len(parts) >= 2:
            taxa.append(parts[0])
            seqs[parts[0]] = parts[-1]
    return taxa, seqs


def _load_ages() -> dict[str, float]:
    df = pd.concat([pd.read_csv(AGES_CASE12), pd.read_csv(AGES_CASE19)], axis=0)
    return {str(r["Name"]).replace("-", "_"): float(r["AGE_SAMPLING"]) for _, r in df.iterrows()}


def _bin_to_nuc(s: str) -> str:
    # Proxy mapping to enable stable BEAST 1.10 run when binary parser backend is unavailable.
    return "".join("A" if c == "0" else ("C" if c == "1" else "N") for c in s)


def _build_xml(taxa: list[str], seqs: dict[str, str], ages: dict[str, float], chain_length: int, stem: str) -> ET.ElementTree:
    beast = ET.Element("beast")

    txs = ET.SubElement(beast, "taxa", {"id": "taxa"})
    for t in taxa:
        tx = ET.SubElement(txs, "taxon", {"id": t})
        ET.SubElement(tx, "date", {"value": f"{ages[t]:.6f}", "direction": "forwards", "units": "years"})

    aln = ET.SubElement(beast, "alignment", {"id": "alignment", "dataType": "nucleotide"})
    for t in taxa:
        s = ET.SubElement(aln, "sequence")
        ET.SubElement(s, "taxon", {"idref": t})
        s.text = f"\n            {_bin_to_nuc(seqs[t])}\n        "

    pats = ET.SubElement(beast, "patterns", {"id": "patterns", "from": "1", "unique": "false"})
    ET.SubElement(pats, "alignment", {"idref": "alignment"})

    init_demo = ET.SubElement(beast, "constantSize", {"id": "initialDemo", "units": "years"})
    ET.SubElement(ET.SubElement(init_demo, "populationSize"), "parameter", {"id": "initialDemo.popSize", "value": "1.0"})
    st = ET.SubElement(beast, "coalescentTree", {"id": "startingTree"})
    ET.SubElement(st, "taxa", {"idref": "taxa"})
    ET.SubElement(st, "constantSize", {"idref": "initialDemo"})

    tm = ET.SubElement(beast, "treeModel", {"id": "treeModel"})
    ET.SubElement(tm, "coalescentTree", {"idref": "startingTree"})
    ET.SubElement(ET.SubElement(tm, "rootHeight"), "parameter", {"id": "treeModel.rootHeight"})
    ET.SubElement(ET.SubElement(tm, "nodeHeights", {"internalNodes": "true"}), "parameter", {"id": "treeModel.internalNodeHeights"})
    ET.SubElement(ET.SubElement(tm, "nodeHeights", {"internalNodes": "true", "rootNode": "true"}), "parameter", {"id": "treeModel.allInternalNodeHeights"})

    hky = ET.SubElement(beast, "hkyModel", {"id": "hky"})
    freq = ET.SubElement(hky, "frequencies")
    fm = ET.SubElement(freq, "frequencyModel", {"dataType": "nucleotide"})
    ET.SubElement(ET.SubElement(fm, "frequencies"), "parameter", {"id": "hky.frequencies", "value": "0.49 0.01 0.49 0.01"})
    ET.SubElement(ET.SubElement(hky, "kappa"), "parameter", {"id": "hky.kappa", "value": "2.0", "lower": "0.0"})

    sm = ET.SubElement(beast, "siteModel", {"id": "siteModel"})
    ET.SubElement(ET.SubElement(sm, "substitutionModel"), "hkyModel", {"idref": "hky"})
    ET.SubElement(ET.SubElement(sm, "mutationRate"), "parameter", {"id": "siteModel.mu", "value": "1.0", "lower": "0.0"})

    br = ET.SubElement(beast, "strictClockBranchRates", {"id": "branchRates"})
    ET.SubElement(ET.SubElement(br, "rate"), "parameter", {"id": "clock.rate", "value": "1.0", "lower": "0.0"})

    tl = ET.SubElement(beast, "treeLikelihood", {"id": "treeLikelihood"})
    ET.SubElement(tl, "patterns", {"idref": "patterns"})
    ET.SubElement(tl, "treeModel", {"idref": "treeModel"})
    ET.SubElement(tl, "siteModel", {"idref": "siteModel"})
    ET.SubElement(tl, "strictClockBranchRates", {"idref": "branchRates"})

    cst = ET.SubElement(beast, "constantSize", {"id": "constant", "units": "years"})
    ET.SubElement(ET.SubElement(cst, "populationSize"), "parameter", {"id": "constant.popSize", "value": "1.0", "lower": "1.0E-8"})
    cl = ET.SubElement(beast, "coalescentLikelihood", {"id": "coalescent"})
    ET.SubElement(ET.SubElement(cl, "model"), "constantSize", {"idref": "constant"})
    ET.SubElement(ET.SubElement(cl, "populationTree"), "treeModel", {"idref": "treeModel"})

    ops = ET.SubElement(beast, "operators", {"id": "operators"})
    for pid in ("hky.kappa", "clock.rate", "treeModel.rootHeight", "constant.popSize"):
        op = ET.SubElement(ops, "scaleOperator", {"scaleFactor": "0.5", "weight": "1"})
        ET.SubElement(op, "parameter", {"idref": pid})
    ET.SubElement(ET.SubElement(ops, "uniformOperator", {"weight": "10"}), "parameter", {"idref": "treeModel.internalNodeHeights"})
    for tag, w in (("subtreeSlide", "5"), ("narrowExchange", "1"), ("wideExchange", "1"), ("wilsonBalding", "1")):
        op = ET.SubElement(ops, tag, {"weight": w} if tag != "subtreeSlide" else {"weight": w, "gaussian": "true", "size": "1.0"})
        ET.SubElement(op, "treeModel", {"idref": "treeModel"})
    de = ET.SubElement(ops, "deltaExchange", {"delta": "0.05", "weight": "1", "autoOptimize": "true"})
    ET.SubElement(de, "parameter", {"idref": "hky.frequencies"})

    mcmc = ET.SubElement(beast, "mcmc", {"id": "mcmc", "chainLength": str(chain_length), "operatorAnalysis": f"{stem}.ops"})
    post = ET.SubElement(mcmc, "posterior", {"id": "posterior"})
    prior = ET.SubElement(post, "prior", {"id": "prior"})
    ET.SubElement(ET.SubElement(prior, "oneOnXPrior"), "parameter", {"idref": "constant.popSize"})
    lik = ET.SubElement(post, "likelihood", {"id": "likelihood"})
    ET.SubElement(lik, "treeLikelihood", {"idref": "treeLikelihood"})
    ET.SubElement(lik, "coalescentLikelihood", {"idref": "coalescent"})
    ET.SubElement(mcmc, "operators", {"idref": "operators"})

    slog = ET.SubElement(mcmc, "log", {"id": "screenLog", "logEvery": "1000"})
    for lbl, tag, ref in [("Posterior", "posterior", "posterior"), ("Prior", "prior", "prior"), ("Likelihood", "likelihood", "likelihood")]:
        c = ET.SubElement(slog, "column", {"label": lbl, "dp": "4", "width": "12"})
        ET.SubElement(c, tag, {"idref": ref})

    flog = ET.SubElement(mcmc, "log", {"id": "fileLog", "logEvery": "1000", "fileName": f"{stem}.log", "overwrite": "true"})
    for tag, ref in [
        ("posterior", "posterior"),
        ("prior", "prior"),
        ("likelihood", "likelihood"),
        ("parameter", "hky.kappa"),
        ("parameter", "clock.rate"),
        ("parameter", "treeModel.rootHeight"),
        ("parameter", "constant.popSize"),
        ("parameter", "hky.frequencies"),
    ]:
        ET.SubElement(flog, tag, {"idref": ref})

    tlog = ET.SubElement(mcmc, "logTree", {"id": "treeFileLog", "logEvery": "1000", "nexusFormat": "true", "fileName": f"{stem}.trees", "sortTranslationTable": "true"})
    ET.SubElement(tlog, "treeModel", {"idref": "treeModel"})
    ET.SubElement(tlog, "posterior", {"idref": "posterior"})

    return ET.ElementTree(beast)


def main() -> None:
    ap = argparse.ArgumentParser(description="Stable BEAST 1.10 run on SNV matrix (HKY proxy mapping 0/1 -> A/C).")
    ap.add_argument("--chain-length", type=int, default=20000)
    ap.add_argument("--stem", default="snv_case12_case19_beast10_hky_proxy")
    args = ap.parse_args()

    taxa, seqs = _parse_nexus(NEXUS_IN)
    ages = _load_ages()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    xml = OUT_DIR / f"{args.stem}.xml"
    _build_xml(taxa, seqs, ages, args.chain_length, args.stem).write(xml, encoding="utf-8", xml_declaration=True)
    cmd = [str(_find_beast10_cmd()), "-overwrite", "-java", str(xml)]
    print("Running:", " ".join(cmd))
    rc = subprocess.call(cmd, cwd=str(OUT_DIR))
    if rc != 0:
        raise SystemExit(rc)
    print(f"Done. Outputs: {OUT_DIR / (args.stem + '.log')} and {OUT_DIR / (args.stem + '.trees')}")


if __name__ == "__main__":
    main()

