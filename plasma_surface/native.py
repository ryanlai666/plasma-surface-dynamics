"""ctypes adapter to an optional locally compiled C++17 backend."""
import ctypes as ct
from functools import lru_cache
import os
from pathlib import Path
import sys
import numpy as np
from .model import Parameters, _check_run, rates


def library_path():
    suffix = ".dll" if os.name == "nt" else ".dylib" if sys.platform == "darwin" else ".so"
    return Path(__file__).resolve().parents[1]/"build"/("surface"+suffix)


@lru_cache(maxsize=1)
def _library():
    path = library_path()
    if not path.is_file():
        raise RuntimeError("C++ backend missing. Run: python scripts/build_cpp.py")
    lib = ct.CDLL(str(path))
    if lib.surface_abi_version() != 1:
        raise RuntimeError("Unsupported C++ ABI; rebuild the backend")
    pointer = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
    base = [pointer, ct.c_int, ct.c_int, ct.c_double]
    lib.surface_mean_field.argtypes = base+[pointer]
    lib.surface_mean_field.restype = ct.c_int
    lib.surface_kmc.argtypes = base+[ct.c_int, ct.c_uint64, ct.c_uint64, pointer]
    lib.surface_kmc.restype = ct.c_int
    return lib


def simulate(phases, p=None, cycles=5, *, method="kmc", sites=256, seed=42, max_events=2_000_000):
    p = p or Parameters()
    _check_run(phases, cycles)
    if method not in ("kmc", "mean_field"):
        raise ValueError("method must be kmc or mean_field")
    if not isinstance(sites, int) or not 1 <= sites <= 2**31-1:
        raise ValueError("sites must be a positive signed 32-bit integer")
    if cycles > 2**31-1 or len(phases) > 2**31-1:
        raise ValueError("Too many cycles or phases for native ABI")
    if not isinstance(seed, int) or not 0 <= seed < 2**64:
        raise ValueError("seed must be an unsigned 64-bit integer")
    if not isinstance(max_events, int) or not 1 <= max_events < 2**53:
        raise ValueError("max_events must be an integer in [1, 2^53)")
    inputs = np.ascontiguousarray([[phase.duration_s, *rates(phase, p)] for phase in phases], dtype=np.float64)
    outputs = np.zeros((cycles*len(phases), 7), dtype=np.float64)
    lib = _library()
    args = [inputs, len(phases), cycles, p.layer_nm]
    if method == "kmc":
        status = lib.surface_kmc(*args, sites, seed, max_events, outputs)
    else:
        status = lib.surface_mean_field(*args, outputs)
    if status:
        raise RuntimeError({1:"Invalid native simulation inputs", 2:"Event budget exceeded", 3:"Native allocation or runtime failure"}.get(status,"Unknown native error"))
    rows = []
    for i, values in enumerate(outputs):
        row = dict(cycle=i//len(phases)+1, phase=phases[i%len(phases)].name,
                   **dict(zip(["time_s", "coverage", "removed_nm", "deposited_nm", "net_removed_nm"], map(float, values[:5]))))
        if method == "kmc":
            row.update(roughness_nm=float(values[5]), events=int(values[6]))
        rows.append(row)
    return rows


def kmc(phases, p=None, cycles=5, sites=256, seed=42, max_events=2_000_000):
    return simulate(phases, p, cycles, sites=sites, seed=seed, max_events=max_events)


def mean_field(phases, p=None, cycles=5):
    return simulate(phases, p, cycles, method="mean_field")
