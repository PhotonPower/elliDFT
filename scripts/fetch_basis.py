#!/usr/bin/env python3
"""Download the basis set files of elliDFT from the Basis Set Exchange.

Writes one NWChem-format file per basis set to data/basis/ and a manifest
(data/basis/MANIFEST.json) with source URL, download date and SHA-256 of each
file. Uses only the Python standard library.

Elements: H to Ar (Z = 1..18). The project scope is H to Cl; Ar is included
because the validation table of the specification lists it as a test atom.

Usage:  python3 scripts/fetch_basis.py [--force]
"""

import argparse
import datetime
import hashlib
import json
import pathlib
import sys
import urllib.parse
import urllib.request

BSE_API = "https://www.basissetexchange.org/api/basis/{name}/format/nwchem/?elements={elements}"
ELEMENTS = list(range(1, 19))

# Basis sets fixed in the specification, section 2. Key: name as used in
# elliDFT, value: role.
BASIS_SETS = {
    "STO-3G": "orbital",
    "6-31G*": "orbital",
    "def2-SVP": "orbital",
    "def2-TZVP": "orbital",
    "def2-SVPD": "orbital",
    "def2-TZVPD": "orbital",
    "cc-pVDZ": "orbital",
    "aug-cc-pVDZ": "orbital",
    "aug-cc-pVTZ": "orbital",
    "def2-universal-JKFIT": "auxiliary",
    "def2-universal-JFIT": "auxiliary",
}


def file_name(basis_name: str) -> str:
    """File name convention shared with src/basis: lower case, '*' -> '_st'."""
    return basis_name.lower().replace("*", "_st") + ".nw"


def fetch(basis_name: str) -> tuple[str, str]:
    url = BSE_API.format(
        name=urllib.parse.quote(basis_name.lower(), safe=""),
        elements=",".join(str(z) for z in ELEMENTS),
    )
    request = urllib.request.Request(url, headers={"User-Agent": "basis-fetch-script"})
    with urllib.request.urlopen(request, timeout=60) as response:
        if response.status != 200:
            raise RuntimeError(f"{url}: HTTP {response.status}")
        return url, response.read().decode("utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--force", action="store_true", help="re-download existing files")
    args = parser.parse_args()

    out_dir = pathlib.Path(__file__).resolve().parent.parent / "data" / "basis"
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "MANIFEST.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}

    for name, role in BASIS_SETS.items():
        path = out_dir / file_name(name)
        if path.exists() and name in manifest and not args.force:
            print(f"keep   {path.name}")
            continue
        url, text = fetch(name)
        if "BASIS" not in text or "END" not in text:
            raise RuntimeError(f"{name}: response does not look like an NWChem basis file")
        path.write_text(text)
        manifest[name] = {
            "file": path.name,
            "role": role,
            "format": "nwchem",
            "elements": "H-Ar (Z=1..18)",
            "source": url,
            "downloaded_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
        }
        print(f"fetch  {path.name}  ({len(text)} bytes)")

    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
