"""Build the optional C++ backend with a local GCC/Clang compiler; no downloads."""

import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
compiler = (
    os.environ.get("CXX") or shutil.which("g++") or shutil.which("clang++") or shutil.which("c++")
)
if not compiler:
    sys.exit("Install a C++17 GCC/Clang compiler or use the Python backend. Set CXX to its path.")
directory = ROOT / "build"
directory.mkdir(exist_ok=True)
suffix = ".dll" if os.name == "nt" else ".dylib" if sys.platform == "darwin" else ".so"
target = directory / ("species" + suffix)
command = [compiler, "-std=c++17", "-O3", "-Wall", "-Wextra", "-shared"]
if os.name == "nt":
    command += ["-static", "-static-libgcc", "-static-libstdc++"]
else:
    command += ["-fPIC"]
command += [str(ROOT / "cpp/species.cpp"), "-o", str(target)]
subprocess.run(command, check=True)
info = dict(
    command=command,
    compiler=subprocess.check_output([compiler, "--version"], text=True).splitlines()[0],
    platform=platform.platform(),
    source_sha256=hashlib.sha256((ROOT / "cpp/species.cpp").read_bytes()).hexdigest(),
    binary_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
)
(directory / "species_build_manifest.json").write_text(
    json.dumps(info, indent=2) + "\n", encoding="utf-8"
)
print(target)
