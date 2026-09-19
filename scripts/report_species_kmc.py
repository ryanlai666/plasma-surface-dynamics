from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
out=ROOT/'docs/species_kmc_results';d=json.loads((out/'summary.json').read_text())
lines=['## Computed conditional values','','| T (K) | Si released / initial motif | SiF4 | SiH2F2 | SiHF3 | NH3 |','|---|---:|---:|---:|---:|---:|']
for r in d['temperatures']:
 g=r['products_per_initial_motif'];lines.append(f"| {r['temperature_K']:.0f} | {r['Si_released_fraction']:.4f} | {g['SiF4']:.4f} | {g['SiH2F2']:.4f} | {g['SiHF3']:.4f} | {g['NH3']:.4f} |")
lines+=['','## Numerical verification','', '- All sampled Si/N/H/F balances and active-center counts are conserved exactly.',f"- Maximum ensemble discrepancy from the master equation: {max(x['max_standard_errors'] for x in d['verification']):.2f} estimated standard errors; acceptance threshold 6.",f"- Independent Python/C++ ensembles: maximum {d['python_cpp_max_standard_errors']:.2f} estimated standard errors in the checked populated states.",'','Five-trial CPU timings below are a local microbenchmark. Random streams differ; event throughput accounts for differing event counts. Concurrent atomistic workloads may affect timings.','','| Backend | Median seconds / 1,000-motif run | Events per second |','|---|---:|---:|']
for name,r in d['benchmark'].items():lines.append(f"| {name} | {r['median_seconds']:.6f} | {r['events_per_second']:.0f} |")
lines+=['','[Full verification, seeds and source hashes](summary.json) | [Compiler manifest](cpp_build.json).','']
p=out/'REPORT.md';p.write_text(p.read_text(encoding='utf-8-sig').split('<!-- RESULTS -->')[0]+'<!-- RESULTS -->\n\n'+'\n'.join(lines),encoding='utf8')
