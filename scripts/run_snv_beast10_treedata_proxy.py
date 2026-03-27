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
    txt = path.read_text(encoding="utf-8", errors="replace")
    m = re.search(r"(?is)\bmatrix\b(.*?)\s*;", txt)
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


def _map_binary_to_nuc(s: str) -> str:
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
        s.text = f"\n            {_map_binary_to_nuc(seqs[t])}\n        "

    pats = ET.SubElement(beast, "patterns", {"id": "patterns", "from": "1", "unique": "false"})
    ET.SubElement(pats, "alignment", {"idref": "alignment"})

    init = ET.SubElement(beast, "constantSize", {"id": "initialDemo", "units": "years"})
    ET.SubElement(ET.SubElement(init, "populationSize"), "parameter", {"id": "initialDemo.popSize", "value": "1.0"})
    start = ET.SubElement(beast, "coalescentTree", {"id": "startingTree"})
    ET.SubElement(start, "taxa", {"idref": "taxa"})
    ET.SubElement(start, "constantSize", {"idref": "initialDemo"})

    tm = ET.SubElement(beast, "treeModel", {"id": "treeModel"})
    ET.SubElement(tm, "coalescentTree", {"idref": "startingTree"})
    ET.SubElement(ET.SubElement(tm, "rootHeight"), "parameter", {"id": "treeModel.rootHeight"})
    ET.SubElement(ET.SubElement(tm, "nodeHeights", {"internalNodes": "true"}), "parameter", {"id": "treeModel.internalNodeHeights"})
    ET.SubElement(ET.SubElement(tm, "nodeHeights", {"internalNodes": "true", "rootNode": "true"}), "parameter", {"id": "treeModel.allInternalNodeHeights"})

    coal = ET.SubElement(beast, "constantSize", {"id": "constant", "units": "years"})
    ET.SubElement(ET.SubElement(coal, "populationSize"), "parameter", {"id": "constant.popSize", "value": "1.0", "lower": "1.0E-8"})
    cl = ET.SubElement(beast, "coalescentLikelihood", {"id": "coalescent"})
    ET.SubElement(ET.SubElement(cl, "model"), "constantSize", {"idref": "constant"})
    ET.SubElement(ET.SubElement(cl, "populationTree"), "treeModel", {"idref": "treeModel"})

    br = ET.SubElement(beast, "strictClockBranchRates", {"id": "branchRates"})
    ET.SubElement(ET.SubElement(br, "rate"), "parameter", {"id": "clock.rate", "value": "1.0", "lower": "0.0"})

    hky = ET.SubElement(beast, "HKYModel", {"id": "hky"})
    fq = ET.SubElement(ET.SubElement(hky, "frequencies"), "frequencyModel", {"dataType": "nucleotide"})
    ET.SubElement(ET.SubElement(fq, "frequencies"), "parameter", {"id": "frequencies", "value": "0.49 0.01 0.49 0.01"})
    ET.SubElement(ET.SubElement(hky, "kappa"), "parameter", {"id": "kappa", "value": "2.0", "lower": "0.0"})

    sm = ET.SubElement(beast, "siteModel", {"id": "siteModel"})
    ET.SubElement(ET.SubElement(sm, "substitutionModel"), "HKYModel", {"idref": "hky"})
    ET.SubElement(ET.SubElement(sm, "relativeRate"), "parameter", {"id": "mu", "value": "1.0", "lower": "0.0"})

    tdl = ET.SubElement(beast, "treeDataLikelihood", {"id": "treeLikelihood", "useAmbiguities": "true"})
    ET.SubElement(tdl, "patterns", {"idref": "patterns"})
    ET.SubElement(tdl, "treeModel", {"idref": "treeModel"})
    ET.SubElement(tdl, "siteModel", {"idref": "siteModel"})
    ET.SubElement(tdl, "strictClockBranchRates", {"idref": "branchRates"})

    ops = ET.SubElement(beast, "operators", {"id": "operators"})
    for pid in ("kappa", "clock.rate", "treeModel.rootHeight", "constant.popSize"):
        op = ET.SubElement(ops, "scaleOperator", {"scaleFactor": "0.5", "weight": "1"})
        ET.SubElement(op, "parameter", {"idref": pid})
    ET.SubElement(ET.SubElement(ops, "uniformOperator", {"weight": "10"}), "parameter", {"idref": "treeModel.internalNodeHeights"})
    for tag in ("subtreeSlide", "narrowExchange", "wideExchange", "wilsonBalding"):
        attrs = {"weight": "5", "gaussian": "true", "size": "1.0"} if tag == "subtreeSlide" else {"weight": "1"}
        op = ET.SubElement(ops, tag, attrs)
        ET.SubElement(op, "treeModel", {"idref": "treeModel"})

    mcmc = ET.SubElement(beast, "mcmc", {"id": "mcmc", "chainLength": str(chain_length), "operatorAnalysis": f"{stem}.ops"})
    post = ET.SubElement(mcmc, "posterior", {"id": "posterior"})
    prior = ET.SubElement(post, "prior", {"id": "prior"})
    ET.SubElement(ET.SubElement(prior, "oneOnXPrior"), "parameter", {"idref": "constant.popSize"})
    lik = ET.SubElement(post, "likelihood", {"id": "likelihood"})
    ET.SubElement(lik, "treeDataLikelihood", {"idref": "treeLikelihood"})
    ET.SubElement(lik, "coalescentLikelihood", {"idref": "coalescent"})
    ET.SubElement(mcmc, "operators", {"idref": "operators"})

    slog = ET.SubElement(mcmc, "log", {"id": "screenLog", "logEvery": "1000"})
    for label, tag, ref in [("Posterior", "posterior", "posterior"), ("Prior", "prior", "prior"), ("Likelihood", "likelihood", "likelihood")]:
        c = ET.SubElement(slog, "column", {"label": label, "dp": "4", "width": "12"})
        ET.SubElement(c, tag, {"idref": ref})

    flog = ET.SubElement(mcmc, "log", {"id": "fileLog", "logEvery": "1000", "fileName": f"{stem}.log", "overwrite": "true"})
    for tag, ref in [
        ("posterior", "posterior"),
        ("prior", "prior"),
        ("likelihood", "likelihood"),
        ("parameter", "kappa"),
        ("parameter", "clock.rate"),
        ("parameter", "treeModel.rootHeight"),
        ("parameter", "constant.popSize"),
    ]:
        ET.SubElement(flog, tag, {"idref": ref})

    tlog = ET.SubElement(mcmc, "logTree", {"id": "treeFileLog", "logEvery": "1000", "nexusFormat": "true", "fileName": f"{stem}.trees", "sortTranslationTable": "true"})
    ET.SubElement(tlog, "treeModel", {"idref": "treeModel"})
    ET.SubElement(tlog, "posterior", {"idref": "posterior"})

    return ET.ElementTree(beast)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain-length", type=int, default=20000)
    ap.add_argument("--stem", default="snv_case12_case19_beast10_treedata_proxy")
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
    print("Done.")


if __name__ == "__main__":
    main()

