import argparse
import json
from pathlib import Path
from .datasets import fetch_hcl, inventory
from .io import manifest, write_csv, write_json
from .workflows import calibrate, campaign, demo


def main():
    parser = argparse.ArgumentParser(description="Plasma Surface Lab: illustrative research workflows")
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("demo", help="Generate ROM, kMC, ML evaluation, and a figure")
    run.add_argument("--output", default="outputs/demo")
    run.add_argument("--samples", type=int, default=128)
    run.add_argument("--seed", type=int, default=42)
    run.add_argument("--backend", choices=["python", "cpp"], default="python")
    sweep = sub.add_parser("sweep", help="Deterministic shardable synthetic ROM campaign")
    sweep.add_argument("--samples", type=int, default=1024)
    sweep.add_argument("--seed", type=int, default=42)
    sweep.add_argument("--shard", type=int, default=0)
    sweep.add_argument("--shards", type=int, default=1)
    sweep.add_argument("--output", default="outputs/sweep.csv")
    sweep.add_argument("--backend", choices=["python", "cpp"], default="python")
    bench = sub.add_parser("benchmark", help="Compare Python and C++ kMC wall times")
    bench.add_argument("--output", default="outputs/benchmark.json")
    bench.add_argument("--sites", type=int, default=2048)
    bench.add_argument("--cycles", type=int, default=10)
    bench.add_argument("--repeats", type=int, default=3)
    materials = sub.add_parser("materials", help="Compare explicit Si/SiNx hypothetical parameter cards")
    materials.add_argument("--cards", default="configs/materials.json")
    materials.add_argument("--output", default="outputs/materials")
    materials.add_argument("--backend", choices=["python", "cpp"], default="python")
    fetch = sub.add_parser("fetch-data", help="Download the 75 kB CC-BY-4.0 Si-HCl DFT archive")
    fetch.add_argument("--output", default="data/raw/si_hcl")
    inspect = sub.add_parser("inventory", help="Inventory a previously downloaded DFT archive")
    inspect.add_argument("archive")
    inspect.add_argument("--output", default="data/processed/structure_inventory.csv")
    fit = sub.add_parser("calibrate", help="Fit two rates to compatible experimental CSV measurements")
    fit.add_argument("csv")
    fit.add_argument("--output", default="outputs/calibration.json")
    args = parser.parse_args()
    if args.command == "demo":
        result = demo(args.output, args.samples, args.seed, args.backend)
    elif args.command == "sweep":
        rows = campaign(args.samples, args.seed, args.shard, args.shards, args.backend)
        if not rows:
            parser.error("Shard has no samples; use fewer shards")
        write_csv(args.output, rows)
        write_json(Path(args.output).with_suffix(".manifest.json"), manifest(vars(args)))
        result = dict(rows=len(rows), output=args.output)
    elif args.command == "benchmark":
        from .benchmark import benchmark
        result = benchmark(args.output, args.sites, args.cycles, args.repeats)
    elif args.command == "materials":
        from .materials import compare
        result = compare(args.cards, args.output, args.backend)
    elif args.command == "fetch-data":
        result = fetch_hcl(args.output)
    elif args.command == "inventory":
        result = inventory(args.archive, args.output)
    else:
        result = calibrate(args.csv, args.output)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
