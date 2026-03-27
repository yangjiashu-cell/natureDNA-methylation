from __future__ import annotations

import html
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DF_ROOT = ROOT / "Duran-FerrerM-evoflux"

HTML_FILES = [
    DF_ROOT / "docs" / "Control_SNPs.html",
    DF_ROOT / "docs" / "SNPs_vs_fCpGs.html",
    DF_ROOT / "docs" / "CNAs_plots.html",
    DF_ROOT / "docs" / "Data_source_Fig.4AB.html",
    DF_ROOT / "docs" / "Nanopore.html",
    DF_ROOT / "docs" / "Nanopore_haplotypes.html",
]

SEARCH_DIRS = [
    ROOT / "CalumGabbutt-evoflux" / "data",
    ROOT / "data",
    ROOT / "evoflux",
    DF_ROOT / "data",
]

REPORT_PATH = ROOT / "figures" / "reproduced" / "step4_genetic_confounding_preflight.md"


def _extract_refs(html_path: Path) -> list[str]:
    txt = html_path.read_text(encoding="utf-8", errors="ignore")
    refs: set[str] = set()
    # Normal quoted paths
    refs.update(
        re.findall(r'"([^"\n<>]*\.(?:csv|tsv|txt|xlsx|RData|gz))"', txt, flags=re.IGNORECASE)
    )
    # HTML escaped quotes
    refs.update(
        re.findall(r"&quot;([^&]+?\.(?:csv|tsv|txt|xlsx|RData|gz))&quot;", txt, flags=re.IGNORECASE)
    )
    # Keep only relative/project-like references (ignore data:image blobs etc)
    keep = []
    for r in refs:
        rr = html.unescape(r).strip()
        if rr.startswith("../") or rr.startswith("../../"):
            keep.append(rr)
    return sorted(set(keep))


def _search_by_name(filename: str) -> list[Path]:
    out: list[Path] = []
    for d in SEARCH_DIRS:
        if not d.exists():
            continue
        out.extend([p for p in d.rglob("*") if p.is_file() and p.name == filename])
    return sorted(set(out))


def main() -> None:
    rows: list[tuple[str, str, str, str]] = []
    total_refs = 0
    direct_ok = 0
    alt_ok = 0

    for hp in HTML_FILES:
        refs = _extract_refs(hp)
        total_refs += len(refs)
        for ref in refs:
            direct = (hp.parent / ref).resolve()
            if direct.exists():
                direct_ok += 1
                rows.append((hp.name, ref, "direct_ok", str(direct)))
                continue

            candidates = _search_by_name(Path(ref).name)
            if candidates:
                alt_ok += 1
                rows.append((hp.name, ref, "alt_name_match", "; ".join(str(c) for c in candidates[:3])))
            else:
                rows.append((hp.name, ref, "missing", "-"))

    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with REPORT_PATH.open("w", encoding="utf-8") as f:
        f.write("# Step4 Genetic Confounding Preflight\n\n")
        f.write("## Scope\n")
        f.write("- HTML workflows checked: Control_SNPs, SNPs_vs_fCpGs, CNAs_plots, Data_source_Fig.4AB, Nanopore, Nanopore_haplotypes.\n")
        f.write("- Data search roots: `CalumGabbutt-evoflux/data`, `data`, `evoflux`, `Duran-FerrerM-evoflux/data`.\n\n")

        f.write("## Summary\n")
        f.write(f"- Total referenced external data paths: **{total_refs}**\n")
        f.write(f"- Directly resolvable from HTML relative paths: **{direct_ok}**\n")
        f.write(f"- Missing direct paths but filename exists somewhere in allowed roots: **{alt_ok}**\n")
        f.write(f"- Unresolved missing dependencies: **{total_refs - direct_ok - alt_ok}**\n\n")

        f.write("## Dependency Matrix\n")
        f.write("| HTML | Referenced path | Status | Resolved candidate |\n")
        f.write("|---|---|---|---|\n")
        for html_name, ref, status, resolved in rows:
            f.write(f"| `{html_name}` | `{ref}` | `{status}` | `{resolved}` |\n")

        f.write("\n## Execution Readiness\n")
        if total_refs - direct_ok - alt_ok == 0:
            f.write("- All referenced data are present by path or safe filename match.\n")
        else:
            f.write("- End-to-end execution from provided HTML code is **not** currently possible due to unresolved inputs.\n")
        f.write("- Also requires an R runtime with packages (`data.table`, `openxlsx`, `minfi`, annotation packages).\n")


if __name__ == "__main__":
    main()

