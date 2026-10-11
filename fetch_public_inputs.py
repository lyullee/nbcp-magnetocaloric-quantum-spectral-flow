"""Download, verify, and optionally extract public research inputs.

Large third-party inputs are intentionally excluded from the repository. This
script uses the frozen manifest and never replaces a file whose hash disagrees.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import tarfile
import urllib.request
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCES = ROOT / "sources"
MANIFEST = SOURCES / "input_manifest_20261011.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def verified(path: Path, expected: str) -> bool:
    return path.is_file() and sha256(path).lower() == expected.lower()


def download(url: str, target: Path, expected: str) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if verified(target, expected):
        print(f"verified existing: {target.relative_to(ROOT)}")
        return
    if target.exists():
        raise RuntimeError(
            f"Existing file has a different hash and was not replaced: {target}"
        )
    temporary = target.with_suffix(target.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as out:
        shutil.copyfileobj(response, out)
    actual = sha256(temporary)
    if actual.lower() != expected.lower():
        temporary.unlink(missing_ok=True)
        raise RuntimeError(f"Hash mismatch for {url}: expected {expected}, got {actual}")
    temporary.replace(target)
    print(f"downloaded and verified: {target.relative_to(ROOT)}")


def safe_extract_tar(archive: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive) as handle:
        handle.extractall(destination, filter="data")
    print(f"extracted: {archive.relative_to(ROOT)} -> {destination.relative_to(ROOT)}")


def safe_extract_zip(archive: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    base = destination.resolve()
    with zipfile.ZipFile(archive) as handle:
        for member in handle.infolist():
            candidate = (destination / member.filename).resolve()
            if candidate != base and base not in candidate.parents:
                raise RuntimeError(f"Unsafe zip member: {member.filename}")
        handle.extractall(destination)
    print(f"extracted: {archive.relative_to(ROOT)} -> {destination.relative_to(ROOT)}")


def clone_at_commit(url: str, destination: Path, commit: str) -> None:
    if destination.exists():
        current = subprocess.check_output(
            ["git", "-C", str(destination), "rev-parse", "HEAD"], text=True
        ).strip()
        if current != commit:
            raise RuntimeError(f"Existing checkout {destination} is at {current}, expected {commit}")
        print(f"verified git checkout: {destination.relative_to(ROOT)} @ {commit}")
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["git", "clone", "--filter=blob:none", url, str(destination)], check=True)
    subprocess.run(["git", "-C", str(destination), "checkout", commit], check=True)
    print(f"cloned: {destination.relative_to(ROOT)} @ {commit}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--download", action="store_true", help="Download missing archive files")
    parser.add_argument("--extract", action="store_true", help="Extract verified tar/zip archives")
    parser.add_argument("--with-qmagen", action="store_true", help="Clone QMagen repositories at frozen commits")
    args = parser.parse_args()

    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    archive_entries = [
        item for item in manifest["files"]
        if item.get("sha256") and not item["file"].endswith("Fig4.pdf")
    ]
    failures = []
    for item in archive_entries:
        target = SOURCES / item["file"]
        if args.download:
            download(item["url"], target, item["sha256"])
        elif verified(target, item["sha256"]):
            print(f"verified existing: {target.relative_to(ROOT)}")
        else:
            failures.append(str(target.relative_to(ROOT)))
            continue

        if args.extract:
            if target.suffix == ".tar":
                safe_extract_tar(target, target.with_suffix(""))
            elif target.suffix == ".zip":
                safe_extract_zip(target, target.with_suffix(""))

    if args.with_qmagen:
        for item in manifest["files"]:
            if "git_commit" in item:
                clone_at_commit(item["url"], SOURCES / item["file"], item["git_commit"])

    if failures:
        print("Missing or hash-mismatched inputs:")
        for name in failures:
            print(f"  {name}")
        print("Run with --download to obtain missing files from their frozen public URLs.")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
