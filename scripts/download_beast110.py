"""
Download official BEAST v1.10.4 (Windows .zip) from GitHub releases into tools/.

Source: https://github.com/beast-dev/beast-mcmc/releases/tag/v1.10.4
Asset: BEAST.v1.10.4.zip

After extraction you should have e.g. tools/BEAST v1.10.4/bin/beast.cmd
(adjust folder name if the zip root differs).

Does not overwrite an existing tools/BEAST v1.10.4 unless --force.
"""
from __future__ import annotations

import argparse
import io
import sys
import zipfile
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
URL = "https://github.com/beast-dev/beast-mcmc/releases/download/v1.10.4/BEAST.v1.10.4.zip"
TARGET_PARENT = ROOT / "tools"
TARGET_DIR = TARGET_PARENT / "BEAST v1.10.4"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true", help="Remove existing TARGET_DIR first")
    args = ap.parse_args()

    if TARGET_DIR.exists() and any(TARGET_DIR.iterdir()):
        if not args.force:
            print(f"Already present (use --force to re-download): {TARGET_DIR}", file=sys.stderr)
            sys.exit(0)
        import shutil

        shutil.rmtree(TARGET_DIR)

    TARGET_PARENT.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {URL} ...")
    with urlopen(URL, timeout=600) as resp:
        data = resp.read()
    zf = zipfile.ZipFile(io.BytesIO(data))
    # Zip often contains a single top-level folder
    zf.extractall(path=TARGET_PARENT)
    zf.close()
    print(f"Extracted under {TARGET_PARENT}")
    # Normalize: if extracted as BEAST v1.10.4 or BEASTv1.10.4, list contents
    for p in TARGET_PARENT.iterdir():
        if "1.10.4" in p.name and p.is_dir():
            beast_cmd = list(p.rglob("beast.cmd"))
            if beast_cmd:
                print(f"Found: {beast_cmd[0]}")
            else:
                print("Warning: beast.cmd not found under extracted tree; check folder layout.", file=sys.stderr)


if __name__ == "__main__":
    main()
