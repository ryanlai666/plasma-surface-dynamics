from pathlib import Path
import platform
import statistics
import time
from . import model, native
from .io import manifest, write_json


def benchmark(output="outputs/benchmark.json", sites=2048, cycles=10, repeats=3):
    if repeats < 1:
        raise ValueError("repeats must be positive")
    phases = model.ale_recipe()
    # Load DLL and exercise both paths before measuring complete Python-call latency.
    model.kmc(phases, cycles=1, sites=16)
    native.kmc(phases, cycles=1, sites=16)
    report = dict(sites=sites, cycles=cycles, repeats=repeats, platform=platform.platform(),
                  note="Same rates and geometry; independent RNG streams. Wall times include Python wrapping. No HPC or GPU measurement.")
    for name, function in [("python", model.kmc), ("cpp", native.kmc)]:
        times, counts, values = [], [], []
        for seed in range(repeats):
            start = time.perf_counter()
            result = function(phases, cycles=cycles, sites=sites, seed=seed)
            times.append(time.perf_counter()-start)
            counts.append(result[-1]["events"])
            values.append(result[-1]["net_removed_nm"])
        report[name] = dict(seconds=times, median_seconds=statistics.median(times), events=counts,
                            net_removed_nm=values, median_events_per_second=statistics.median([n/t for n,t in zip(counts,times)]))
    report["wall_time_speedup"] = report["python"]["median_seconds"]/report["cpp"]["median_seconds"]
    report["provenance"] = manifest(dict(sites=sites, cycles=cycles, repeats=repeats))
    import json
    build = native.library_path().parent/"build_manifest.json"
    report["native_build"] = json.loads(build.read_text()) if build.exists() else None
    write_json(output, report)
    return report
