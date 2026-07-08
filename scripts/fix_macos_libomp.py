#!/usr/bin/env python3
"""Repoint duplicate bundled libomp.dylib copies (torch, scikit-learn, ...) to a
single canonical libomp.dylib so only one OpenMP runtime is ever loaded into the
process. This is what actually causes the

    OMP: Error #15: Initializing libomp.dylib, but found libomp.dylib already
    initialized.

crash on macOS -- several wheels each bundle their own copy of libomp.dylib, and
loading more than one into the same process trips libomp's safety check. Rather
than silencing that check with KMP_DUPLICATE_LIB_OK=TRUE, this makes sure there is
only one libomp.dylib to load in the first place.

Run once after installing boltz (pip install -e .):

    python scripts/fix_macos_libomp.py

Safe to re-run; already-canonical binaries are left untouched.
"""

import subprocess
import sys
import sysconfig
from pathlib import Path


def find_canonical(site_packages: Path) -> Path:
    conda_prefix = subprocess.os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        candidate = Path(conda_prefix) / "lib" / "libomp.dylib"
        if candidate.exists():
            return candidate
    candidate = site_packages / "torch" / "lib" / "libomp.dylib"
    if candidate.exists():
        return candidate
    msg = "Could not find a libomp.dylib to use as the canonical copy."
    raise SystemExit(msg)


def otool_deps(path: Path) -> list[str]:
    out = subprocess.run(
        ["otool", "-L", str(path)], capture_output=True, text=True, check=True
    ).stdout
    # Skip the first line: it's the library's own install name (or the path, for
    # executables), not a dependency.
    return [line.split()[0] for line in out.splitlines()[1:]]


def resolve(dep: str, loader_dir: Path) -> Path | None:
    if dep.startswith("@loader_path/"):
        return (loader_dir / dep.removeprefix("@loader_path/")).resolve()
    if dep.startswith("@rpath/"):
        # Already resolved relative to the binary's rpath entries (typically the
        # canonical env lib dir already) -- nothing to fix here.
        return None
    if dep.startswith("/"):
        return Path(dep).resolve()
    return None


def main() -> None:
    site_packages = Path(sysconfig.get_paths()["purelib"])
    canonical = find_canonical(site_packages).resolve()
    print(f"Canonical libomp.dylib: {canonical}")

    binaries = [
        p
        for p in site_packages.rglob("*")
        if p.suffix in (".dylib", ".so") and p.resolve() != canonical
    ]

    changed = 0
    for binary in binaries:
        if binary.name == "libomp.dylib":
            continue  # don't rewrite other libomp copies' own install name
        try:
            deps = otool_deps(binary)
        except subprocess.CalledProcessError:
            continue

        for dep in deps:
            if Path(dep).name != "libomp.dylib":
                continue
            resolved = resolve(dep, binary.parent)
            if resolved is None or not resolved.exists():
                continue
            if resolved == canonical:
                continue

            print(f"Repointing {binary}\n  {dep} -> {canonical}")
            subprocess.run(
                ["install_name_tool", "-change", dep, str(canonical), str(binary)],
                check=True,
            )
            # install_name_tool invalidates the code signature; macOS refuses to
            # load unsigned/mis-signed binaries, so re-sign ad hoc.
            subprocess.run(
                ["codesign", "--force", "-s", "-", str(binary)],
                check=True,
                capture_output=True,
            )
            changed += 1

    print(f"Done. Repointed {changed} binaries.")


if __name__ == "__main__":
    if sys.platform != "darwin":
        msg = "This script is only needed on macOS."
        raise SystemExit(msg)
    main()
