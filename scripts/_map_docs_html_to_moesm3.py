from __future__ import annotations

import html
import json
import re
from pathlib import Path


def main() -> None:
    docs_dir = Path(r"d:\naturedna\evoflux-reproduce\Duran-FerrerM-evoflux\docs")
    moesm3_path = Path(r"d:\naturedna\evoflux-reproduce\evoflux\41586_2025_9374_MOESM3_ESM.xlsx")

    # Patterns for explicit file loaders inside the knitted HTML R code blocks.
    patterns = [
        re.compile(r'read\.xlsx\(\s*(?:xlsxFile\s*=\s*)?"([^"]+)"'),
        re.compile(r'openxlsx::read\.xlsx\(\s*"([^"]+)"'),
        re.compile(r'fread\(\s*"([^"]+)"'),
        re.compile(r'read\.csv\(\s*"([^"]+)"'),
        re.compile(r'load\(\s*"([^"]+)"'),
        re.compile(r'readRDS\(\s*"([^"]+)"'),
        re.compile(r'getSheetNames\(\s*"([^"]+)"'),
    ]

    html_files = sorted(docs_dir.glob("*.html"))
    if not html_files:
        raise SystemExit(f"No HTML files found in {docs_dir}")

    mapping: dict[str, dict[str, object]] = {}
    for f in html_files:
        txt = html.unescape(f.read_text(encoding="utf-8", errors="ignore"))
        hits: list[str] = []
        for pat in patterns:
            hits.extend(m.group(1) for m in pat.finditer(txt))
        # Preserve order, de-dup
        uniq: list[str] = []
        seen = set()
        for h in hits:
            if h not in seen:
                seen.add(h)
                uniq.append(h)

        # Heuristic: does the HTML mention supplementary_table_* (MOESM3-like)?
        mentions_moesm3_style = bool(
            re.search(r"supplementary_table_\d+", txt, flags=re.IGNORECASE)
        )

        mapping[f.name] = {
            "explicit_loaded_files": uniq,
            "mentions_supplementary_table_x": mentions_moesm3_style,
        }

    out = {
        "docs_dir": str(docs_dir),
        "moesm3_path": str(moesm3_path),
        "html_count": len(html_files),
        "html_to_inputs": mapping,
    }

    print(json.dumps(out, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

