"""
Fig.4 PISCA (BEAST 1.8.4 + PISCA) pipeline using Zenodo harmonised fCpG matrix.

Inputs (under repo root):
  - data/beta_fcpgs.csv       — 978 fCpGs × 2204 samples (rows = CpG IDs)
  - data/QC_2204_metadata.csv — AGE_SAMPLING and sample IDs

Discrete states (per sample):
  - gmm: 3-component Gaussian mixture on beta (sorted ascending mean → 0,1,2).
  - beta_em: 3-component Beta mixture via EM (scipy); closer to paper than GMM,
    still not identical to full Bayesian Stan.
  - stan: optional; requires `cmdstanpy` + CmdStan + `stan/beta_mixture_3state.stan`,
    or run `scripts/fit_beta_mixture_stan.py` and pass `--discretize external`.
  - external: CSV matrix (samples × CpG IDs) with states 0/1/2 from your own Stan run.
  - threshold: β<0.2→0, β>0.8→2 else 1 (fast fallback).

Population / path sampling:
  - Default coalescent: constant population size (matches PISCA template).
  - `--population-model exponential`: exponential growth coalescent (for model comparison
    with path sampling / marginal likelihood in BEAST 1.8.4; see `marginalLikelihoodEstimator`
    in examples). `--write-path-sampling-appendix` writes a fragment README under
    `figures/reproduced/pisca_zenodo_fig4/` describing how to run step sampling (not executed here).

Priors (closer to Nature 2025 Methods than PISCA example XML):
  - Relative substitution rates: LogNormal mean=1, stdev=0.6, meanInRealSpace=true
  - Clock rate: Exponential(mean≈0.104) as a positive prior in the same spirit as
    HalfNormal(0, 0.13) (exact half-normal not available as a single BEAST element).

MCMC:
  - Default 2_000_000 generations, logEvery 2000 (edit via CLI).
  - Paper: 1e8 generations, logEvery 1e5 — use --paper-chain (very long run).

Paper-aligned shortcut:
  - `--paper-methods` sets: `--discretize stan`, `--paper-chain`, and writes path-sampling notes
    (same Stan chain length / BEAST log as paper-scale; very long runtime).
  - Stan per-sample fit uses the same warmup/sampling defaults as `scripts/fit_beta_mixture_stan.py`
    (500 warmup, 500 sampling) unless you override via code constants below.

Outputs:
  figures/reproduced/pisca_zenodo_fig4/<stem>.xml|.log|.trees|.ops
  + TreeAnnotator MCC: <stem>_MCC.tree (10% burn-in, median heights)
"""
from __future__ import annotations

import argparse
import math
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BETA_PATH = ROOT / "data" / "beta_fcpgs.csv"
META_PATH = ROOT / "data" / "QC_2204_metadata.csv"
BEAST_CMD = ROOT / "tools" / "BEAST v1.8.4" / "bin" / "beast.cmd"
TREEANNOTATOR_CMD = ROOT / "tools" / "BEAST v1.8.4" / "bin" / "treeannotator.cmd"
TEMPLATE_XML = ROOT / "tools" / "PISCAv1.1" / "examples" / "biallelicBinary4Params.xml"
OUT_DIR = ROOT / "figures" / "reproduced" / "pisca_zenodo_fig4"

# Match `scripts/fit_beta_mixture_stan.py` (paper Methods: Bayesian Beta-mixture in Stan).
STAN_ITER_WARMUP = 500
STAN_ITER_SAMPLING = 500
STAN_CHAINS = 2

# Longitudinal tips (Zenodo IDs). Case 12: T1–T4 without duplicate T1 aliases.
CASE12_SAMPLES = ["SCLL-545", "SCLL-546", "SCLL-547", "SCLL-548"]
# Case 19: T1,T2,T3,T4,T5 (SCLL-018 omitted as near-duplicate of SCLL-531).
CASE19_SAMPLES = ["SCLL-493", "SCLL-494", "SCLL-531", "SCLL-532", "SCLL-533"]


def _beta_to_state_threshold(v: float) -> str:
    if v < 0.2:
        return "0"
    if v > 0.8:
        return "2"
    return "1"


def _discretize_sample_gmm(beta: np.ndarray, random_state: int = 0) -> str:
    from sklearn.mixture import GaussianMixture

    x = np.asarray(beta, dtype=float).reshape(-1, 1)
    x = np.clip(x, 1e-4, 1.0 - 1e-4)
    gm = GaussianMixture(
        n_components=3,
        covariance_type="diag",
        reg_covar=1e-4,
        random_state=random_state,
        max_iter=500,
    )
    gm.fit(x)
    means = gm.means_.ravel()
    order = np.argsort(means)
    inv = np.empty_like(order)
    inv[order] = np.arange(3)
    labels = gm.predict(x)
    states = inv[labels]
    return "".join(str(int(s)) for s in states)


def _discretize_sample_beta_em(beta: np.ndarray, random_state: int = 0) -> str:
    """3-component Beta mixture EM; order components by ascending E[Beta]=alpha/(alpha+beta)."""
    from scipy.optimize import minimize
    from scipy.stats import beta as beta_dist

    rng = np.random.RandomState(random_state)
    x = np.clip(np.asarray(beta, dtype=float).ravel(), 1e-4, 1.0 - 1e-4)
    n = x.size
    k_comp = 3

    def weighted_beta_nll(par: np.ndarray, xs: np.ndarray, w: np.ndarray) -> float:
        a, b = float(par[0]), float(par[1])
        if a <= 0.05 or b <= 0.05:
            return 1e200
        ll = beta_dist.logpdf(xs, a, b)
        return float(-np.sum(w * ll))

    def fit_component(xs: np.ndarray, w: np.ndarray) -> tuple[float, float]:
        m = np.average(xs, weights=w)
        v = np.average((xs - m) ** 2, weights=w)
        v = max(v, 1e-8)
        t = m * (1 - m) / v - 1
        if t <= 0:
            a0, b0 = 1.0, 1.0
        else:
            a0 = max(0.1, m * t)
            b0 = max(0.1, (1 - m) * t)
        res = minimize(
            weighted_beta_nll,
            np.array([a0, b0], dtype=float),
            args=(xs, w),
            method="L-BFGS-B",
            bounds=((0.05, 80.0), (0.05, 80.0)),
        )
        return float(res.x[0]), float(res.x[1])

    # Init responsibilities with random soft assignment
    z = rng.dirichlet(np.ones(k_comp), size=n)
    pi = z.mean(axis=0)
    alpha = np.ones(k_comp)
    beta_p = np.ones(k_comp)

    for _ in range(150):
        # M-step
        for k in range(k_comp):
            alpha[k], beta_p[k] = fit_component(x, z[:, k])
        pi = z.sum(axis=0) / n
        pi = np.clip(pi, 1e-8, None)
        pi /= pi.sum()
        # E-step
        log_rho = np.zeros((n, k_comp))
        for k in range(k_comp):
            log_rho[:, k] = np.log(pi[k]) + beta_dist.logpdf(x, alpha[k], beta_p[k])
        log_rho -= log_rho.max(axis=1, keepdims=True)
        rho = np.exp(log_rho)
        rho /= rho.sum(axis=1, keepdims=True)
        z = rho

    means_ord = alpha / (alpha + beta_p)
    order = np.argsort(means_ord)
    inv = np.empty(k_comp, dtype=int)
    inv[order] = np.arange(k_comp)
    states = np.argmax(z, axis=1)
    states = inv[states]
    return "".join(str(int(s)) for s in states)


def _load_external_state_matrix(
    path: Path,
    samples: list[str],
    cpg_ids: list[str],
) -> pd.DataFrame:
    ext = pd.read_csv(path, index_col=0)
    ext.index = ext.index.astype(str)
    ext.columns = ext.columns.astype(str)
    cpg_ids = [str(x) for x in cpg_ids]
    missing_s = set(samples) - set(ext.index)
    if missing_s:
        raise ValueError(f"external CSV missing samples: {sorted(missing_s)}")
    missing_c = set(cpg_ids) - set(ext.columns)
    if missing_c:
        raise ValueError(
            f"external CSV missing {len(missing_c)} CpG columns (example): {sorted(missing_c)[:5]}"
        )
    sub = ext.loc[samples, cpg_ids].astype(int)
    if sub.min().min() < 0 or sub.max().max() > 2:
        raise ValueError("external states must be integers 0, 1, or 2")
    return sub


def _discretize_sample_stan_cmdstan(beta: np.ndarray, stan_file: Path) -> str:
    try:
        import cmdstanpy
    except ImportError as e:
        raise RuntimeError(
            "discretize=stan requires cmdstanpy and CmdStan (set CMDSTAN). "
            "Or use scripts/fit_beta_mixture_stan.py and --discretize external."
        ) from e
    if not stan_file.exists():
        raise FileNotFoundError(stan_file)
    y = np.clip(np.asarray(beta, dtype=float).ravel(), 1e-4, 1.0 - 1e-4)
    sm = cmdstanpy.CmdStanModel(stan_file=str(stan_file))
    fit = sm.sample(
        data={"N": int(y.size), "y": y.tolist()},
        chains=STAN_CHAINS,
        iter_warmup=STAN_ITER_WARMUP,
        iter_sampling=STAN_ITER_SAMPLING,
        show_progress=False,
        show_console=False,
    )
    from scipy.stats import beta as beta_dist

    draws = fit.draws_pd()
    means: list[float] = []
    for k in (1, 2, 3):
        a = float(draws[f"alpha[{k}]"].median())
        b = float(draws[f"beta[{k}]"].median())
        means.append(a / (a + b))
    order = np.argsort(means)
    inv = np.empty(3, dtype=int)
    inv[order] = np.arange(3)
    ta = [float(draws[f"alpha[{k}]"].median()) for k in (1, 2, 3)]
    tb = [float(draws[f"beta[{k}]"].median()) for k in (1, 2, 3)]
    tt = [float(draws[f"theta[{k}]"].median()) for k in (1, 2, 3)]
    s = sum(tt)
    tt = [t / s for t in tt]
    chars: list[str] = []
    for val in y:
        logp = [np.log(tt[k]) + beta_dist.logpdf(val, ta[k], tb[k]) for k in range(3)]
        k_hat = int(np.argmax(logp))
        chars.append(str(int(inv[k_hat])))
    return "".join(chars)


def _load_beta_and_meta() -> tuple[pd.DataFrame, pd.DataFrame]:
    if not BETA_PATH.exists():
        raise FileNotFoundError(f"Missing Zenodo beta matrix: {BETA_PATH}")
    if not META_PATH.exists():
        raise FileNotFoundError(f"Missing metadata: {META_PATH}")
    beta = pd.read_csv(BETA_PATH, index_col=0)
    if beta.shape[0] != 978:
        raise ValueError(f"Expected 978 fCpG rows, got {beta.shape[0]}")
    meta = pd.read_csv(META_PATH)
    return beta, meta


def _meta_for_samples(samples: list[str], meta: pd.DataFrame) -> pd.DataFrame:
    sub = meta[meta["PARTICIPANT_ID_ANONYMOUS"].isin(samples)][
        ["PARTICIPANT_ID_ANONYMOUS", "AGE_SAMPLING"]
    ].drop_duplicates(subset=["PARTICIPANT_ID_ANONYMOUS"])
    missing = set(samples) - set(sub["PARTICIPANT_ID_ANONYMOUS"])
    if missing:
        raise ValueError(f"Samples missing from QC_2204_metadata.csv: {sorted(missing)}")
    sub = sub.set_index("PARTICIPANT_ID_ANONYMOUS").loc[samples].reset_index()
    sub = sub.rename(columns={"PARTICIPANT_ID_ANONYMOUS": "Name"})
    return sub


def _build_state_matrix(
    samples: list[str],
    beta: pd.DataFrame,
    discretize: str,
    *,
    external_csv: Path | None = None,
    stan_file: Path | None = None,
    random_state: int = 0,
) -> pd.DataFrame:
    cols = [c for c in samples if c in beta.columns]
    if len(cols) != len(samples):
        missing = set(samples) - set(beta.columns)
        raise ValueError(f"Samples missing from beta_fcpgs.csv: {sorted(missing)}")
    cpg_ids = list(beta.index.astype(str))
    if discretize == "external":
        if not external_csv or not external_csv.exists():
            raise ValueError("discretize=external requires --external-states-csv pointing to a CSV")
        return _load_external_state_matrix(external_csv, samples, cpg_ids).astype(str)

    sub = beta[cols].T.astype(float)
    rows: list[list[str]] = []
    for sid in samples:
        b = sub.loc[sid].values
        if discretize == "threshold":
            rows.append([_beta_to_state_threshold(float(v)) for v in b])
        elif discretize == "gmm":
            s = _discretize_sample_gmm(b, random_state=random_state)
            rows.append(list(s))
        elif discretize == "beta_em":
            s = _discretize_sample_beta_em(b, random_state=random_state)
            rows.append(list(s))
        elif discretize == "stan":
            sf = stan_file or (ROOT / "stan" / "beta_mixture_3state.stan")
            s = _discretize_sample_stan_cmdstan(b, sf)
            rows.append(list(s))
        else:
            raise ValueError(discretize)
    mat = pd.DataFrame(rows, index=samples, columns=beta.index)
    return mat


def _replace_taxa_and_alignment(
    root: ET.Element,
    state_mat: pd.DataFrame,
    ages: pd.DataFrame,
) -> None:
    taxa = root.find("taxa")
    alignment = root.find("alignment")
    if taxa is None or alignment is None:
        raise RuntimeError("Template missing <taxa> or <alignment>")

    for child in list(taxa):
        taxa.remove(child)
    for child in list(alignment):
        alignment.remove(child)

    ET.SubElement(alignment, "dataType").set("idref", "biallelicBinary")

    age_map = dict(zip(ages["Name"], ages["AGE_SAMPLING"].astype(float)))

    for sample in state_mat.index:
        age = float(age_map[sample])
        tid = sample.replace("-", "_")
        taxon = ET.SubElement(taxa, "taxon")
        taxon.set("id", tid)
        date = ET.SubElement(taxon, "date")
        date.set("value", f"{age:.6f}")
        date.set("direction", "forwards")
        date.set("units", "years")

        seq = "".join(state_mat.loc[sample].values)
        sequence = ET.SubElement(alignment, "sequence")
        taxon_ref = ET.SubElement(sequence, "taxon")
        taxon_ref.set("idref", tid)
        taxon_ref.tail = f"\n\t\t{seq}\n\t"


def _strip_operators_for_small_trees(operators: ET.Element, n_taxa: int) -> None:
    if n_taxa >= 3:
        return
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


def _apply_population_model(root: ET.Element, mode: str) -> None:
    """Swap constant-size coalescent for exponential growth (BEAST 1.8.4 XML)."""
    if mode == "constant":
        return
    if mode != "exponential":
        raise ValueError(mode)

    beast = root
    const_el = None
    const_idx = None
    for i, el in enumerate(beast):
        if el.tag == "constantSize" and el.get("id") == "constant":
            const_el = el
            const_idx = i
            break
    if const_el is None or const_idx is None:
        raise RuntimeError("Template missing <constantSize id='constant'>")

    beast.remove(const_el)
    exp = ET.Element("exponentialGrowth")
    exp.set("id", "exponential")
    exp.set("units", "years")
    ps = ET.SubElement(exp, "populationSize")
    p0 = ET.SubElement(ps, "parameter")
    p0.set("id", "exponential.popSize")
    p0.set("value", "1")
    p0.set("lower", "0.0")
    gr = ET.SubElement(exp, "growthRate")
    p1 = ET.SubElement(gr, "parameter")
    p1.set("id", "exponential.growthRate")
    p1.set("value", "0.0")
    p1.set("lower", "-Infinity")
    p1.set("upper", "Infinity")
    beast.insert(const_idx, exp)

    for sim in root.findall(".//coalescentSimulator"):
        for ch in list(sim):
            if ch.tag == "constantSize" and ch.get("idref") == "constant":
                sim.remove(ch)
                ref = ET.Element("exponentialGrowth")
                ref.set("idref", "exponential")
                sim.append(ref)

    for cl in root.findall(".//coalescentLikelihood"):
        model = cl.find("model")
        if model is None:
            continue
        for ch in list(model):
            if ch.tag == "constantSize" and ch.get("idref") == "constant":
                model.remove(ch)
                ref = ET.Element("exponentialGrowth")
                ref.set("idref", "exponential")
                model.append(ref)

    operators = root.find("operators")
    if operators is not None:
        replaced = False
        for ch in list(operators):
            if ch.tag == "scaleOperator" and ch.find("parameter[@idref='constant.popSize']") is not None:
                operators.remove(ch)
                s1 = ET.Element("scaleOperator")
                s1.set("scaleFactor", "0.5")
                s1.set("weight", "3.0")
                p = ET.Element("parameter")
                p.set("idref", "exponential.popSize")
                s1.append(p)
                operators.append(s1)
                rw = ET.Element("randomWalkOperator")
                rw.set("windowSize", "1.0")
                rw.set("weight", "3.0")
                p2 = ET.Element("parameter")
                p2.set("idref", "exponential.growthRate")
                rw.append(p2)
                operators.append(rw)
                replaced = True
                break
        if not replaced:
            raise RuntimeError("No scaleOperator on constant.popSize to swap for exponential model")

    prior = root.find(".//prior[@id='prior']")
    if prior is None:
        raise RuntimeError("No <prior id='prior'>")
    for el in list(prior):
        if el.tag == "oneOnXPrior":
            pr = el.find("parameter[@idref='constant.popSize']")
            if pr is not None:
                pr.set("idref", "exponential.popSize")
    lap = ET.Element("laplacePrior")
    lap.set("mean", "0.0")
    lap.set("scale", "10.0")
    lp = ET.Element("parameter")
    lp.set("idref", "exponential.growthRate")
    lap.append(lp)
    prior.append(lap)

    for fl in root.findall(".//log[@id='fileLog']"):
        has_growth = False
        for pr in fl.findall("parameter"):
            if pr.get("idref") == "constant.popSize":
                pr.set("idref", "exponential.popSize")
            if pr.get("idref") == "exponential.growthRate":
                has_growth = True
        if not has_growth:
            g = ET.Element("parameter")
            g.set("idref", "exponential.growthRate")
            fl.append(g)


def _write_path_sampling_appendix(out_dir: Path) -> Path:
    """Text-only note; marginal likelihood is not run automatically."""
    p = out_dir / "path_sampling_beast18_notes.txt"
    text = """Path sampling / stepping-stone (BEAST 1.8.4)

Paper Methods: compare demographic models (e.g. constant vs exponential population size)
using marginal likelihood / path sampling.

Steps (manual; not executed by this repo):
1. Generate two XMLs with this script: --population-model constant and
   --population-model exponential (same discretization, priors, chain length).
2. Add a marginalLikelihoodEstimator block to each XML (see BEAST 1.8.4 example
   tools/BEAST v1.8.4/examples/testXML/testTreePathSampling.xml).
3. Run path sampling / stepping-stone as in that example; compare log marginal
   likelihoods in Tracer or from BEAST stdout.

- PISCA uses a cenancestor likelihood; verify compatibility with your BEAST build
  before relying on marginal likelihood estimates.
- This repository does not execute path sampling automatically (runtime + plugin risk).
"""
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")
    return p


def _apply_paper_priors_and_heights(root: ET.Element, max_tip_age: float) -> None:
    """Set luca_height / luca_branch bounds; swap priors for clock + relative rates."""
    luca_h = max_tip_age + 5.0
    for p in root.findall(".//parameter"):
        pid = p.get("id")
        if pid == "luca_height":
            p.set("value", f"{luca_h:.6f}")
        elif pid == "luca_branch":
            p.set("value", "5")
            p.set("upper", "250")
            p.set("lower", "0.0")
        elif pid == "clock.rate":
            p.set("value", "0.02")

    for u in root.findall(".//uniformPrior"):
        if u.find(".//parameter[@idref='luca_branch']") is not None:
            u.set("upper", "250")
            u.set("lower", "0.0")

    prior = root.find(".//prior[@id='prior']")
    if prior is None:
        raise RuntimeError("No <prior id='prior'>")

    to_remove: list[ET.Element] = []
    for el in list(prior):
        if el.tag == "logNormalPrior":
            ref = el.find("parameter[@idref='clock.rate']")
            if ref is not None:
                to_remove.append(el)
        if el.tag == "exponentialPrior":
            ch = list(el)
            if ch and ch[0].get("idref") in {
                "biallelicBinary.demethylation",
                "biallelicBinary.homozygousMethylation",
                "biallelicBinary.homozygousDemethylation",
            }:
                to_remove.append(el)
    for el in to_remove:
        prior.remove(el)

    # Half-normal(0, 0.13): mean = σ*sqrt(2/pi) ≈ 0.1039 — exponential with this mean.
    exp_mean = 0.13 * math.sqrt(2.0 / math.pi)
    exp_clock = ET.Element("exponentialPrior")
    exp_clock.set("mean", f"{exp_mean:.6f}")
    exp_clock.set("offset", "0.0")
    exp_clock.append(_prior_param_ref("clock.rate"))

    # Insert clock prior after oneOnXPrior (before relative-rate priors)
    insert_at = next(i for i, c in enumerate(prior) if c.tag == "oneOnXPrior") + 1
    prior.insert(insert_at, exp_clock)

    for pid in (
        "biallelicBinary.demethylation",
        "biallelicBinary.homozygousMethylation",
        "biallelicBinary.homozygousDemethylation",
    ):
        ln = ET.Element("logNormalPrior")
        ln.set("mean", "1")
        ln.set("stdev", "0.6")
        ln.set("offset", "0.0")
        ln.set("meanInRealSpace", "true")
        ln.append(_prior_param_ref(pid))
        prior.append(ln)


def _prior_param_ref(idref: str) -> ET.Element:
    p = ET.Element("parameter")
    p.set("idref", idref)
    return p


def _configure_mcmc(
    root: ET.Element,
    stem: str,
    chain_length: int,
    log_every: int,
    n_taxa: int,
) -> None:
    mcmc = root.find("mcmc")
    if mcmc is None:
        raise RuntimeError("No <mcmc>")
    mcmc.set("chainLength", str(int(chain_length)))
    mcmc.set("operatorAnalysis", f"{stem}.ops")

    operators = root.find("operators")
    if operators is not None:
        _strip_operators_for_small_trees(operators, n_taxa)

    for node in root.findall(".//log"):
        if node.get("id") == "screenLog":
            node.set("logEvery", str(log_every))
        if node.get("id") == "fileLog":
            node.set("logEvery", str(log_every))
            node.set("fileName", f"{stem}.log")
            node.set("overwrite", "true")
    for node in root.findall(".//logTree"):
        node.set("logEvery", str(log_every))
        node.set("fileName", f"{stem}.trees")


def _write_xml(tree: ET.ElementTree, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tree.write(path, encoding="utf-8", xml_declaration=True)


def _count_trees_in_nexus(p: Path) -> int:
    text = p.read_text(encoding="utf-8", errors="replace")
    return sum(1 for line in text.splitlines() if line.strip().startswith("tree "))


def _run_treeannotator_mcc(trees_path: Path, stem: str) -> Path:
    n = _count_trees_in_nexus(trees_path)
    burnin = max(1, int(0.10 * n))
    out = trees_path.with_name(f"{stem}_MCC.tree")
    cmd = [
        str(TREEANNOTATOR_CMD),
        "-burninTrees",
        str(burnin),
        "-heights",
        "median",
        str(trees_path),
        str(out),
    ]
    r = subprocess.run(cmd, cwd=str(trees_path.parent), check=False)
    if r.returncode != 0:
        raise RuntimeError(f"TreeAnnotator failed: {cmd}")
    return out


def run_case(
    stem: str,
    samples: list[str],
    *,
    discretize: str,
    chain_length: int,
    log_every: int,
    run_beast: bool,
    run_annotator: bool,
    population_model: str,
    external_states_csv: Path | None,
    stan_file: Path | None,
    random_state: int,
    also_write_exponential_xml: bool = False,
) -> None:
    beta, meta_all = _load_beta_and_meta()
    ages = _meta_for_samples(samples, meta_all)
    state_mat = _build_state_matrix(
        samples,
        beta,
        discretize,
        external_csv=external_states_csv,
        stan_file=stan_file,
        random_state=random_state,
    )
    max_age = float(ages["AGE_SAMPLING"].max())

    tree = ET.parse(TEMPLATE_XML)
    root = tree.getroot()
    _replace_taxa_and_alignment(root, state_mat, ages)
    _apply_population_model(root, population_model)
    _apply_paper_priors_and_heights(root, max_age)
    _configure_mcmc(root, stem, chain_length, log_every, n_taxa=len(samples))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    xml_path = OUT_DIR / f"{stem}.xml"
    _write_xml(tree, xml_path)

    if also_write_exponential_xml and population_model == "constant":
        tree_e = ET.parse(xml_path)
        root_e = tree_e.getroot()
        _apply_population_model(root_e, "exponential")
        exp_stem = f"{stem}_exponential"
        _configure_mcmc(root_e, exp_stem, chain_length, log_every, n_taxa=len(samples))
        exp_path = OUT_DIR / f"{exp_stem}.xml"
        _write_xml(tree_e, exp_path)
        print(f"Wrote exponential coalescent duplicate for path sampling: {exp_path}")

    tag = discretize if discretize != "external" else "external"
    state_mat.to_csv(OUT_DIR / f"{stem}_states_{tag}.csv")
    ages.to_csv(OUT_DIR / f"{stem}_ages.csv", index=False)

    if run_beast and not BEAST_CMD.exists():
        raise FileNotFoundError(f"BEAST not found: {BEAST_CMD}")

    if run_beast:
        cmd = [str(BEAST_CMD), "-overwrite", "-beagle_off", str(xml_path)]
        print("Running:", " ".join(cmd))
        r = subprocess.run(cmd, cwd=str(OUT_DIR), check=False)
        if r.returncode != 0:
            raise RuntimeError(f"BEAST exited with {r.returncode}")
        if run_annotator:
            tp = OUT_DIR / f"{stem}.trees"
            if not tp.exists():
                raise FileNotFoundError(tp)
            mcc = _run_treeannotator_mcc(tp, stem)
            print(f"MCC tree: {mcc}")


def main() -> None:
    ap = argparse.ArgumentParser(description="PISCA Fig.4 pipeline (Zenodo matrix).")
    ap.add_argument(
        "--cases",
        default="case12,case19",
        help="Comma-separated: case12, case19",
    )
    ap.add_argument(
        "--discretize",
        choices=("gmm", "threshold", "beta_em", "stan", "external"),
        default="gmm",
        help="Per-sample state coding: beta_em (Beta-mix EM, closer to paper than gmm), gmm, threshold, stan, external.",
    )
    ap.add_argument(
        "--external-states-csv",
        type=Path,
        default=None,
        help="With --discretize external: CSV rows=samples, cols=CpG IDs, values 0/1/2.",
    )
    ap.add_argument(
        "--stan-file",
        type=Path,
        default=None,
        help="Override path to beta_mixture_3state.stan (default: repo stan/).",
    )
    ap.add_argument("--random-state", type=int, default=0, help="Seed for gmm / beta_em.")
    ap.add_argument(
        "--population-model",
        choices=("constant", "exponential"),
        default="constant",
        help="Coalescent demographic model (exponential for comparison with path sampling).",
    )
    ap.add_argument(
        "--write-path-sampling-appendix",
        action="store_true",
        help="Write path_sampling_beast18_notes.txt next to outputs (no MCMC run).",
    )
    ap.add_argument("--chain-length", type=int, default=2_000_000)
    ap.add_argument("--log-every", type=int, default=2000)
    ap.add_argument(
        "--paper-chain",
        action="store_true",
        help="Use paper-like 1e8 states / 1e5 sampling (very long).",
    )
    ap.add_argument(
        "--paper-methods",
        action="store_true",
        help="Use Stan discretization + paper chain length + path-sampling notes (Nature 2025 Methods).",
    )
    ap.add_argument(
        "--also-write-exponential-xml",
        action="store_true",
        help=(
            "When --population-model constant: also write <stem>_exponential.xml (exponential growth) "
            "for marginal likelihood / path sampling comparison (same priors & states)."
        ),
    )
    ap.add_argument("--skip-beast", action="store_true", help="Only write XML + CSV previews.")
    ap.add_argument("--no-annotator", action="store_true", help="Skip TreeAnnotator after BEAST.")
    args = ap.parse_args()

    if args.paper_methods:
        args.discretize = "stan"
        args.paper_chain = True
        args.write_path_sampling_appendix = True

    if args.paper_chain:
        args.chain_length = 100_000_000
        args.log_every = 100_000

    if args.discretize == "external" and not args.external_states_csv:
        raise SystemExit("--discretize external requires --external-states-csv")

    case_map = {
        "case12": ("fig4a_case12_zenodo", CASE12_SAMPLES),
        "case19": ("fig4b_case19_zenodo", CASE19_SAMPLES),
    }
    for key in args.cases.split(","):
        key = key.strip()
        if key not in case_map:
            raise SystemExit(f"Unknown case {key!r}; choose from {list(case_map)}")
        stem, samples = case_map[key]
        run_case(
            stem,
            samples,
            discretize=args.discretize,
            chain_length=args.chain_length,
            log_every=args.log_every,
            run_beast=not args.skip_beast,
            run_annotator=not args.skip_beast and not args.no_annotator,
            population_model=args.population_model,
            external_states_csv=args.external_states_csv,
            stan_file=args.stan_file,
            random_state=args.random_state,
            also_write_exponential_xml=args.also_write_exponential_xml,
        )
        print(f"Done {stem}")
    if args.write_path_sampling_appendix:
        p = _write_path_sampling_appendix(OUT_DIR)
        print(f"Wrote {p}")


if __name__ == "__main__":
    main()
