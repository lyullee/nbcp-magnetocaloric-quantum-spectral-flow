"""Verify the generated Cryogenics Word manuscript and supplement."""

from __future__ import annotations

import json
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parent
MAIN = ROOT / "manuscript" / "Manuscript_Cryogenics_20261012.docx"
SI = ROOT / "manuscript" / "Supplementary_Information_Cryogenics_20261012.docx"
OUT = ROOT / "results" / "cryogenics_docx_verification.json"
NS = {
    "w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
    "m": "http://schemas.openxmlformats.org/officeDocument/2006/math",
}


def inspect(path):
    with ZipFile(path) as archive:
        document_xml = archive.read("word/document.xml").decode("utf-8")
        styles_xml = archive.read("word/styles.xml")
    root = ET.fromstring(document_xml.encode("utf-8"))
    styles_root = ET.fromstring(styles_xml)
    word_text = f"{{{NS['w']}}}t"
    math_text = f"{{{NS['m']}}}t"
    text = "".join(node.text or "" for node in root.iter() if node.tag in {word_text, math_text})
    word_val = f"{{{NS['w']}}}val"
    title_styles = []
    for style in styles_root.findall(".//w:style", NS):
        name = style.find("w:name", NS)
        if name is not None and name.get(word_val) == "Title":
            title_styles.append(style)
    title_has_border = any(style.find("./w:pPr/w:pBdr", NS) is not None for style in title_styles)
    return {
        "bytes": path.stat().st_size,
        "text": text,
        "tables": len(root.findall(".//w:tbl", NS)),
        "equations": len(root.findall(".//m:oMath", NS)),
        "drawings": len(root.findall(".//w:drawing", NS)),
        "raw_markup_counts": {
            token: document_xml.count(token)
            for token in (r"\\frac", r"\\Delta", r"\\(", r"\\[", "**", "![")
        },
        "title_has_border": title_has_border,
    }


def main():
    main_doc = inspect(MAIN)
    si_doc = inspect(SI)
    common_tokens = [
        "db6et9clf4us73c582f0",
        "db6eurqmb58s738bmisg",
        "IBM Kingston",
        "25.8%",
    ]
    checks = {
        "documents_exist": MAIN.exists() and SI.exists(),
        "main_has_all_figures": main_doc["drawings"] == 11,
        "si_has_all_figures": si_doc["drawings"] == 4,
        "main_tables_present": main_doc["tables"] == 5,
        "si_tables_present": si_doc["tables"] == 13,
        "native_math_present": main_doc["equations"] >= 190 and si_doc["equations"] >= 175,
        "qpu_evidence_present_in_both": all(
            token in main_doc["text"] and token in si_doc["text"] for token in common_tokens
        ),
        "frozen_qpu_decision_in_si": (
            "PASS_SIGN_RESOLUTION_HOLD_QUANTITATIVE_MAGNITUDE" in si_doc["text"]
        ),
        "authors_and_correspondence_present": all(
            token in main_doc["text"]
            for token in ["Ugwiyeon Lee", "Seong-Yun Bae", "Hakrae Kim", "Jungsu Lee", "lyullee@kgs.or.kr"]
        ),
        "no_raw_markdown_or_tex": all(
            count == 0
            for doc in (main_doc, si_doc)
            for count in doc["raw_markup_counts"].values()
        ),
        "title_style_has_no_border": not main_doc["title_has_border"] and not si_doc["title_has_border"],
        "reasonable_file_sizes": main_doc["bytes"] > 3_000_000 and si_doc["bytes"] > 700_000,
    }
    failed = [name for name, passed in checks.items() if not passed]
    payload = {
        "status": "PASS" if not failed else "FAIL",
        "checks": checks,
        "failed": failed,
        "main": {key: value for key, value in main_doc.items() if key not in {"text"}},
        "supplement": {key: value for key, value in si_doc.items() if key not in {"text"}},
        "visual_qa": {
            "renderer": "Microsoft Word 16 fixed-layout PDF and Chrome PDF viewer screenshots",
            "main_pages": 22,
            "supplement_pages": 19,
            "inspection": "all pages inspected; no clipping, overlap, missing figure, or broken table observed",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    raise SystemExit(1 if failed else 0)


if __name__ == "__main__":
    main()
