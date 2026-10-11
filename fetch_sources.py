"""Download only the public source workbooks linked by Xiang et al. (2024)."""
from pathlib import Path
import concurrent.futures
import hashlib
import json
import urllib.request

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "sources"
RAW.mkdir(exist_ok=True)
BASE = "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-023-06885-w/MediaObjects/"
SOURCES = [
    ("figure_2.xlsx", "41586_2023_6885_MOESM2_ESM.xlsx", "Main figure 2: cooling and heat capacity"),
    ("extended_figure_4.xlsx", "41586_2023_6885_MOESM6_ESM.xlsx", "Extended figure 4: nuclear/phonon corrections and entropy"),
    ("extended_figure_6.xlsx", "41586_2023_6885_MOESM7_ESM.xlsx", "Extended figure 6: refrigeration comparisons and hold times"),
]

def fetch(item):
    local_name, remote_name, description = item
    url = BASE + remote_name
    target = RAW / local_name
    if target.exists():
        data = target.read_bytes()
    else:
        with urllib.request.urlopen(url, timeout=50) as response:
            data = response.read()
        if not data.startswith(b"PK"):
            raise ValueError(f"Not an XLSX ZIP: {local_name}")
        target.write_bytes(data)
    return {"file": local_name, "url": url, "description": description,
            "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}

with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
    entries = list(executor.map(fetch, SOURCES))
manifest = {
    "access_date": "2026-10-11", "paper_doi": "10.1038/s41586-023-06885-w",
    "paper_url": "https://www.nature.com/articles/s41586-023-06885-w",
    "model_source": "https://arxiv.org/html/2504.11298v1#S0.SS5",
    "files": entries,
    "note": "Public source-data files; not a claim of access to every original laboratory log. Source files are preserved unchanged."
}
(ROOT / "source_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(manifest, ensure_ascii=False, indent=2))
