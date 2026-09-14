"""Command-line entry points; scientific inputs are never silently fabricated."""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
from .data import load_dataset, load_csv
from .demo import make_demo
from .pipeline import analyze, update_manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description="Local ET kinetics and DOS/lambda/coupling decomposition")
    commands = parser.add_subparsers(dest="command", required=True)
    convert = commands.add_parser("import-csv", help="validate and join site/LDOS CSV exports")
    convert.add_argument("--sites", type=Path, required=True)
    convert.add_argument("--ldos", type=Path, required=True)
    convert.add_argument("--metadata", type=Path, required=True)
    convert.add_argument("--output", type=Path, required=True)
    demo = commands.add_parser("demo", help="run explicitly synthetic moire examples")
    demo.add_argument("--angles", nargs="+", type=float, default=[0.8, 1.1, 1.5, 2.0])
    demo.add_argument("--grid-size", type=int, default=24)
    demo.add_argument("--energy-points", type=int, default=1601)
    external = commands.add_parser("analyze", help="analyze independently computed arrays")
    external.add_argument("--input", type=Path, required=True)
    external.add_argument("--metadata", type=Path, required=True)
    external.add_argument("--reference", required=True, help="explicit reference site ID")
    external.add_argument("--target", required=True, help="explicit comparison site ID")
    for command in (demo, external):
        command.add_argument("--output", type=Path, required=True, help="new or empty output directory")
        command.add_argument("--temperature", type=float, default=298.15)
        command.add_argument("--eta", type=float, default=0.0)
        command.add_argument("--mu", type=float, default=0.0)
        command.add_argument("--area", type=float, default=1.0, help="effective coupling area in nm²")
        command.add_argument("--direction", choices=["reduction", "oxidation"], default="reduction")
        command.add_argument("--no-plots", action="store_true")
    args = parser.parse_args(argv)
    if args.command == "import-csv":
        try:
            data = load_csv(args.sites, args.ldos, args.metadata)
            if args.output.exists() and any(args.output.iterdir()):
                raise ValueError("output is not empty; choose a new directory")
            data.save(args.output)
            update_manifest(args.output)
            print(f"Validated {len(data.site_ids)} sites: {args.output}")
            return
        except (ValueError, OSError, KeyError, TypeError) as exc:
            parser.exit(2, f"tbg-et: {exc}\n")
    options = dict(temperature_k=args.temperature, eta_v=args.eta, mu_ev=args.mu,
                   area_nm2=args.area, direction=args.direction, plots=not args.no_plots)
    try:
        if args.command == "analyze":
            result = analyze(load_dataset(args.input, args.metadata), args.output,
                             args.reference, args.target, **options)
            print(json.dumps({"output": str(args.output), "input_kind": result["input_kind"],
                              "log10_rate_ratio": result["decomposition"]["log10_rate_ratio"]}))
        else:
            if args.output.exists() and any(args.output.iterdir()):
                raise ValueError("output is not empty; choose a new directory")
            angles = sorted(set(args.angles))
            names = [f"theta_{a:g}" for a in angles]
            if len(set(names)) != len(names):
                raise ValueError("angles are too close to have distinct output names")
            # Validate all angles and common parameters before creating any output.
            datasets = [make_demo(a, args.grid_size, args.energy_points, args.temperature) for a in angles]
            rows = []
            for data, name, theta in zip(datasets, names, angles):
                result = analyze(data, args.output / name, "AB-center", "AA-center", **options)
                factors = result["decomposition"]["shapley_log_factors"]
                rows.append({"theta_deg": theta, **{k: v / np.log(10) for k, v in factors.items()},
                             "total": result["decomposition"]["log10_rate_ratio"]})
            with (args.output / "twist_sweep.csv").open("w", newline="") as file:
                writer = csv.DictWriter(file, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            if not args.no_plots:
                from .plotting import plot_twist
                plot_twist(rows, args.output / "twist_sweep.png")
            update_manifest(args.output)
            print(f"SYNTHETIC demo complete: {args.output}. No experimental or first-principles claim.")
    except (ValueError, OSError, KeyError, TypeError) as exc:
        parser.exit(2, f"tbg-et: {exc}\n")


if __name__ == "__main__":
    main()
