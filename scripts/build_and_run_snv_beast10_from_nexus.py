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
    candidates = [
        ROOT / "tools" / "BEAST v1.10.4" / "bin" / "beast.cmd",
        ROOT / "tools" / "BEAST v1.10.0" / "bin" / "beast.cmd",
        ROOT / "tools" / "BEAST v1.10.5" / "bin" / "beast.cmd",
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError("BEAST 1.10.x beast.cmd not found under tools/")


def _parse_nexus_matrix(path: Path) -> tuple[list[str], dict[str, str]]:
    txt = path.read_text(encoding="utf-8", errors="replace")
    m = re.search(r"(?is)\bmatrix\b(.*?)\s*;", txt)
    if not m:
        raise ValueError(f"No MATRIX block in {path}")
    block = m.group(1).strip()
    taxa: list[str] = []
    seqs: dict[str, str] = {}
    for ln in block.splitlines():
        s = ln.strip()
        if not s:
            continue
        parts = re.split(r"\s+", s)
        if len(parts) < 2:
            continue
        tid, seq = parts[0], parts[-1]
        taxa.append(tid)
        seqs[tid] = seq
    if not taxa:
        raise ValueError("No taxa parsed from MATRIX")
    return taxa, seqs


def _load_ages() -> dict[str, float]:
    a12 = pd.read_csv(AGES_CASE12)
    a19 = pd.read_csv(AGES_CASE19)
    ages = pd.concat([a12, a19], axis=0, ignore_index=True)
    out: dict[str, float] = {}
    for _, r in ages.iterrows():
        # NEXUS taxa use underscore form (SCLL_545)
        out[str(r["Name"]).replace("-", "_")] = float(r["AGE_SAMPLING"])
    return out


def _build_xml(taxa: list[str], seqs: dict[str, str], ages: dict[str, float], chain_length: int, stem: str) -> ET.ElementTree:
    beast = ET.Element("beast")

    taxa_el = ET.SubElement(beast, "taxa", {"id": "taxa"})
    for t in taxa:
        tx = ET.SubElement(taxa_el, "taxon", {"id": t})
        d = ET.SubElement(tx, "date")
        d.set("value", f"{ages[t]:.6f}")
        d.set("direction", "forwards")
        d.set("units", "years")

    gdt = ET.SubElement(beast, "generalDataType", {"id": "binaryDataType"})
    ET.SubElement(gdt, "state", {"code": "0"})
    ET.SubElement(gdt, "state", {"code": "1"})
    ET.SubElement(gdt, "ambiguity", {"code": "?", "states": "01"})
    ET.SubElement(gdt, "ambiguity", {"code": "-", "states": "01"})

    aln = ET.SubElement(beast, "alignment", {"id": "alignment"})
    ET.SubElement(aln, "dataType", {"idref": "binaryDataType"})
    for t in taxa:
        seq = ET.SubElement(aln, "sequence")
        ET.SubElement(seq, "taxon", {"idref": t})
        seq.text = f"\n            {seqs[t]}\n        "

    pats = ET.SubElement(beast, "patterns", {"id": "patterns", "from": "1", "unique": "false"})
    ET.SubElement(pats, "alignment", {"idref": "alignment"})

    init_demo = ET.SubElement(beast, "constantSize", {"id": "initialDemo", "units": "years"})
    ps = ET.SubElement(init_demo, "populationSize")
    ET.SubElement(ps, "parameter", {"id": "initialDemo.popSize", "value": "1.0", "lower": "1.0E-8"})

    stree = ET.SubElement(beast, "coalescentTree", {"id": "startingTree"})
    ET.SubElement(stree, "taxa", {"idref": "taxa"})
    ET.SubElement(stree, "constantSize", {"idref": "initialDemo"})

    tm = ET.SubElement(beast, "treeModel", {"id": "treeModel"})
    ET.SubElement(tm, "coalescentTree", {"idref": "startingTree"})
    rh = ET.SubElement(tm, "rootHeight")
    ET.SubElement(rh, "parameter", {"id": "treeModel.rootHeight"})
    nh1 = ET.SubElement(tm, "nodeHeights", {"internalNodes": "true"})
    ET.SubElement(nh1, "parameter", {"id": "treeModel.internalNodeHeights"})
    nh2 = ET.SubElement(tm, "nodeHeights", {"internalNodes": "true", "rootNode": "true"})
    ET.SubElement(nh2, "parameter", {"id": "treeModel.allInternalNodeHeights"})

    fm = ET.SubElement(beast, "frequencyModel", {"id": "binaryFreqModel"})
    ET.SubElement(fm, "dataType", {"idref": "binaryDataType"})
    fqs = ET.SubElement(fm, "frequencies")
    ET.SubElement(fqs, "parameter", {"id": "binary.frequencies", "dimension": "2", "value": "0.5 0.5"})

    gsm = ET.SubElement(beast, "generalSubstitutionModel", {"id": "binaryModel"})
    sf = ET.SubElement(gsm, "frequencies")
    ET.SubElement(sf, "frequencyModel", {"idref": "binaryFreqModel"})
    rates = ET.SubElement(gsm, "rates", {"relativeTo": "1"})
    ET.SubElement(rates, "parameter", {"id": "binary.rates", "dimension": "0"})

    sm = ET.SubElement(beast, "siteModel", {"id": "siteModel"})
    subm = ET.SubElement(sm, "substitutionModel")
    ET.SubElement(subm, "generalSubstitutionModel", {"idref": "binaryModel"})
    mr = ET.SubElement(sm, "mutationRate")
    ET.SubElement(mr, "parameter", {"id": "siteModel.mu", "value": "1.0", "lower": "0.0"})

    br = ET.SubElement(beast, "strictClockBranchRates", {"id": "branchRates"})
    rate = ET.SubElement(br, "rate")
    ET.SubElement(rate, "parameter", {"id": "clock.rate", "value": "1.0", "lower": "0.0"})

    tl = ET.SubElement(beast, "treeLikelihood", {"id": "treeLikelihood"})
    ET.SubElement(tl, "patterns", {"idref": "patterns"})
    ET.SubElement(tl, "treeModel", {"idref": "treeModel"})
    ET.SubElement(tl, "siteModel", {"idref": "siteModel"})
    ET.SubElement(tl, "strictClockBranchRates", {"idref": "branchRates"})

    cst = ET.SubElement(beast, "constantSize", {"id": "constant", "units": "years"})
    cps = ET.SubElement(cst, "populationSize")
    ET.SubElement(cps, "parameter", {"id": "constant.popSize", "value": "1.0", "lower": "1.0E-8"})

    coal = ET.SubElement(beast, "coalescentLikelihood", {"id": "coalescent"})
    ET.SubElement(coal, "model").append(ET.Element("constantSize", {"idref": "constant"}))
    ET.SubElement(coal, "populationTree").append(ET.Element("treeModel", {"idref": "treeModel"}))

    ops = ET.SubElement(beast, "operators", {"id": "operators"})
    op1 = ET.SubElement(ops, "scaleOperator", {"scaleFactor": "0.5", "weight": "1"})
    ET.SubElement(op1, "parameter", {"idref": "clock.rate"})
    op2 = ET.SubElement(ops, "scaleOperator", {"scaleFactor": "0.5", "weight": "1"})
    ET.SubElement(op2, "parameter", {"idref": "treeModel.rootHeight"})
    op3 = ET.SubElement(ops, "uniformOperator", {"weight": "10"})
    ET.SubElement(op3, "parameter", {"idref": "treeModel.internalNodeHeights"})
    op4 = ET.SubElement(ops, "subtreeSlide", {"weight": "5", "gaussian": "true", "size": "1.0"})
    ET.SubElement(op4, "treeModel", {"idref": "treeModel"})
    op5 = ET.SubElement(ops, "narrowExchange", {"weight": "1"})
    ET.SubElement(op5, "treeModel", {"idref": "treeModel"})
    op6 = ET.SubElement(ops, "wideExchange", {"weight": "1"})
    ET.SubElement(op6, "treeModel", {"idref": "treeModel"})
    op7 = ET.SubElement(ops, "wilsonBalding", {"weight": "1"})
    ET.SubElement(op7, "treeModel", {"idref": "treeModel"})
    op8 = ET.SubElement(ops, "scaleOperator", {"scaleFactor": "0.5", "weight": "2"})
    ET.SubElement(op8, "parameter", {"idref": "constant.popSize"})
    op9 = ET.SubElement(ops, "deltaExchange", {"delta": "0.1", "weight": "1", "autoOptimize": "true"})
    ET.SubElement(op9, "parameter", {"idref": "binary.frequencies"})

    mcmc = ET.SubElement(beast, "mcmc", {"id": "mcmc", "chainLength": str(chain_length), "operatorAnalysis": f"{stem}.ops"})
    post = ET.SubElement(mcmc, "posterior", {"id": "posterior"})
    prior = ET.SubElement(post, "prior", {"id": "prior"})
    ET.SubElement(prior, "oneOnXPrior").append(ET.Element("parameter", {"idref": "constant.popSize"}))
    lik = ET.SubElement(post, "likelihood", {"id": "likelihood"})
    ET.SubElement(lik, "treeLikelihood", {"idref": "treeLikelihood"})
    ET.SubElement(lik, "coalescentLikelihood", {"idref": "coalescent"})

    ET.SubElement(mcmc, "operators", {"idref": "operators"})

    slog = ET.SubElement(mcmc, "log", {"id": "screenLog", "logEvery": "1000"})
    for label, ref in [("Posterior", "posterior"), ("Prior", "prior"), ("Likelihood", "likelihood"), ("RootHeight", "treeModel.rootHeight")]:
        c = ET.SubElement(slog, "column", {"label": label, "dp": "4", "width": "12"})
        tag = "parameter" if ref == "treeModel.rootHeight" else ("posterior" if ref == "posterior" else ("prior" if ref == "prior" else "likelihood"))
        ET.SubElement(c, tag, {"idref": ref})

    flog = ET.SubElement(mcmc, "log", {"id": "fileLog", "logEvery": "1000", "fileName": f"{stem}.log", "overwrite": "true"})
    for ref in ["posterior", "prior", "likelihood", "clock.rate", "treeModel.rootHeight", "constant.popSize", "binary.frequencies"]:
        if ref == "posterior":
            tag = "posterior"
        elif ref == "prior":
            tag = "prior"
        elif ref == "likelihood":
            tag = "likelihood"
        else:
            tag = "parameter"
        ET.SubElement(flog, tag, {"idref": ref})

    tlog = ET.SubElement(mcmc, "logTree", {"id": "treeFileLog", "logEvery": "1000", "nexusFormat": "true", "fileName": f"{stem}.trees", "sortTranslationTable": "true"})
    ET.SubElement(tlog, "treeModel", {"idref": "treeModel"})
    ET.SubElement(tlog, "posterior", {"idref": "posterior"})

    rep = ET.SubElement(beast, "report")
    prop = ET.SubElement(rep, "property", {"name": "timer"})
    ET.SubElement(prop, "mcmc", {"idref": "mcmc"})

    return ET.ElementTree(beast)


def main() -> None:
    ap = argparse.ArgumentParser(description="Build + optionally run BEAST1.10 on existing SNV nexus.")
    ap.add_argument("--nexus", type=Path, default=NEXUS_IN)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--stem", default="snv_case12_case19_beast10")
    ap.add_argument("--chain-length", type=int, default=20000)
    ap.add_argument("--skip-run", action="store_true")
    args = ap.parse_args()

    taxa, seqs = _parse_nexus_matrix(args.nexus)
    ages = _load_ages()
    missing = [t for t in taxa if t not in ages]
    if missing:
        raise ValueError(f"Missing ages for taxa: {missing}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    xml_path = args.out_dir / f"{args.stem}.xml"
    tree = _build_xml(taxa, seqs, ages, args.chain_length, args.stem)
    tree.write(xml_path, encoding="utf-8", xml_declaration=True)
    print(f"Wrote XML: {xml_path}")

    if args.skip_run:
        return

    beast_cmd = _find_beast10_cmd()
    # BEAST 1.10 does not accept -beagle_off; default backend selection is sufficient here.
    cmd = [str(beast_cmd), "-overwrite", "-java", str(xml_path)]
    print("Running:", " ".join(cmd))
    rc = subprocess.call(cmd, cwd=str(args.out_dir))
    if rc != 0:
        raise SystemExit(rc)
    print("Done.")


if __name__ == "__main__":
    main()

