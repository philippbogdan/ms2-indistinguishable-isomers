"""Validate the evidence and export eight CSV tables for the four figures."""

import argparse
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ms2.results import calculate, check_expected


def make_data(source, output, expected):
    source, output = Path(source), Path(output)
    result = calculate(source)
    if not check_expected(result, expected)["matches_expected"]:
        raise ValueError("The evidence differs from expected.json; do not render it")
    units = pd.read_parquet(source / "isomer_units.parquet")
    repeats = pd.read_parquet(source / "repeat_spectra.parquet")
    files = {
        "isomer_scores.csv": units[[
            "query_id", "reference_molecule", "score", "instrument_type", "precursor_mz",
            "fragment_peak_count", "collision_energy", "identical_peak_list", "exact_peak_list",
        ]],
        "repeat_scores.csv": repeats[["query_id", "fold", "score"]],
        "threshold.csv": pd.DataFrame([{
            "quantile": result["method"]["repeat_quantile"], "threshold": result["threshold"],
            "isomer_units": result["units"]["n"], "within_units": result["units"]["within"],
            "within_percent": result["units"]["percent"],
            "training_repeats": result["repeats"]["train"]["n"],
        }]),
    }
    for name in ["example_pairs.csv", "example_spectra.csv", "example_peaks.csv"]:
        files[name] = pd.read_csv(source / name, float_precision="round_trip")
    files["example_spectra.csv"] = files["example_spectra.csv"][[
        "pair_id", "role", "spectrum_id", "name", "molecule_key", "smiles", "precursor_mz",
    ]]
    groups = {
        "instrument": ["Orbitrap", "QTOF"],
        "collision_energy": ["0", "(0,20]", "(20,30]", "(30,40]", "(40,60]", "over 60"],
        "fragment_peaks": ["0 to 2", "3 to 5", "6 to 10", "11 to 20", "21 to 50", "over 50"],
        "precursor_mz": ["up to 200", "(200,300]", "(300,400]", "(400,500]", "over 500"],
    }
    rows = []
    for dimension, labels in groups.items():
        for order, label in enumerate(labels):
            if dimension == "collision_energy" and label == "0":
                values = result["breakdowns"]["collision_energy_outside_bands"]
                if values["n"] != result["collision_energy_zero_units"]:
                    raise ValueError("Unexpected energies outside the positive bands")
            else:
                values = result["breakdowns"][dimension][label]
            rows.append({"dimension": dimension, "order": order, "label": label,
                         **{key: values[key] for key in ["n", "within", "percent"]}})
    files["breakdowns.csv"] = pd.DataFrame(rows)
    files["sensitivity.csv"] = pd.DataFrame([
        {"quantile": float(label.rstrip("%")) / 100,
         **{key: values[key] for key in ["threshold", "n", "within", "percent"]}}
        for label, values in result["sensitivity"].items()
    ])
    output.mkdir(parents=True, exist_ok=True)
    unexpected = {path.name for path in output.iterdir()} - set(files)
    if unexpected:
        raise ValueError(f"Unexpected files in the data-only destination: {sorted(unexpected)}")
    for name, frame in files.items():
        frame.to_csv(output / name, index=False)
    print(f"Wrote {len(files)} checked data files to {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data", help="Repository data directory")
    parser.add_argument("--output", type=Path, default=ROOT / "figures" / "data", help="Directory for the eight generated CSV tables")
    parser.add_argument("--expected", type=Path, default=ROOT / "data" / "expected.json",
                        help="Fixed reference values")
    args = parser.parse_args()
    make_data(args.data, args.output, args.expected)


if __name__ == "__main__":
    main()
