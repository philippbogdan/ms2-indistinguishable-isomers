"""Regenerate every reported count and share from the included data, offline."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .examples import validate_examples
from .protocol import (
    DATA_PATH, DUPLICATE_ATOL, DUPLICATE_RTOL, EXAMPLE_IDS, INPUT_BYTES,
    INTENSITY_POWER, PRECURSOR_MARGIN_DA, QUANTILE, REVISION, ROOT, SHA256,
    TOLERANCE_DA, TOLERANCE_PPM, URL,
)


def share(numerator, denominator):
    """Keep counts with each fraction so its denominator is always explicit."""
    numerator, denominator = int(numerator), int(denominator)
    fraction = numerator / denominator if denominator else None
    return {"n": denominator, "within": numerator, "fraction": fraction,
            "percent": None if fraction is None else 100 * fraction}


def grouped_shares(frame, labels):
    return {
        str(label): share(group.within.sum(), len(group))
        for label, group in frame.groupby(labels, observed=True, sort=True)
    }


def pair_summary(frame, columns):
    pairs = frame.groupby(columns, sort=False).within.agg(["sum", "size"])
    return {
        "n": len(pairs),
        "at_least_one": share((pairs["sum"] > 0).sum(), len(pairs)),
        "at_least_half": share((2 * pairs["sum"] >= pairs["size"]).sum(), len(pairs)),
        "all": share((pairs["sum"] == pairs["size"]).sum(), len(pairs)),
    }


def agreement(frame):
    """Count ordered repeat pairs algebraically, without expanding their rows."""
    groups = frame.groupby(["query_molecule", "reference_molecule", "condition"]).within.agg(["sum", "size"])
    groups = groups[groups["size"] >= 2]
    inside, total = groups["sum"], groups["size"]
    outside = total - inside
    return {
        "given_within": share((inside * (inside - 1)).sum(), (inside * (total - 1)).sum()),
        "given_outside": share((outside * (outside - 1)).sum(), (outside * (total - 1)).sum()),
    }


def validate_tables(units, repeats, spectra, molecules):
    """Fail on broken joins, conditions or folds before reporting a headline."""
    if not spectra.spectrum_id.is_unique:
        raise ValueError("Duplicate spectrum identifiers")
    meta = spectra.set_index("spectrum_id", drop=False)
    if units.duplicated(["query_id", "reference_molecule"]).any():
        raise ValueError("Duplicate isomer units")
    if not repeats.query_id.is_unique:
        raise ValueError("Duplicate repeat queries")
    for frame, isomer in [(units, True), (repeats, False)]:
        query = meta.loc[frame.query_id].reset_index(drop=True)
        reference = meta.loc[frame.reference_id].reset_index(drop=True)
        witness_column = "same_condition_reference_id" if isomer else "same_condition_repeat_id"
        witness = meta.loc[frame[witness_column]].reset_index(drop=True)
        expected_key = frame.reference_molecule if isomer else frame.molecule_key
        checks = [
            (query.adduct == reference.adduct).all(),
            (reference.molecule_key == expected_key).all(),
            (witness.molecule_key == expected_key).all(),
            (query.condition == witness.condition).all(),
            (query.condition == frame.condition).all(),
            query.instrument_type.notna().all(), query.collision_energy.notna().all(),
            np.isfinite(frame.score).all(), (frame.score >= 0).all(),
            (frame.score <= 1 + 1e-12).all(),
            ((query.condition == reference.condition) == frame.reference_same_condition).all(),
        ]
        if isomer:
            checks += [
                (query.molecule_key == frame.query_molecule).all(),
                (query.molecule_key != reference.molecule_key).all(),
                (query.formula == frame.formula).all(), (witness.formula == frame.formula).all(),
                (query.fold == frame.query_fold).all(), (reference.fold == frame.reference_fold).all(),
                (query.instrument_type == frame.instrument_type).all(),
                (query.collision_energy == frame.collision_energy).all(),
                (query.precursor_mz == frame.precursor_mz).all(),
                (query.fragment_peak_count == frame.fragment_peak_count).all(),
                (~frame.exact_peak_list | frame.identical_peak_list).all(),
            ]
        else:
            checks += [(query.molecule_key == frame.molecule_key).all(),
                       (query.spectrum_id != reference.spectrum_id).all(),
                       (query.spectrum_id != witness.spectrum_id).all(),
                       (query.fold == frame.fold).all()]
        if not all(checks):
            raise ValueError(f"Inconsistent {'isomer' if isomer else 'repeat'} evidence table")
    expected_molecules = spectra[["molecule_key", "smiles", "formula", "fold"]].drop_duplicates()
    if set(map(tuple, expected_molecules.to_numpy())) != set(map(tuple, molecules.to_numpy())):
        raise ValueError("Molecule records do not match the spectrum metadata")
    if molecules.groupby("molecule_key").fold.nunique().max() != 1:
        raise ValueError("A molecule occurs in more than one fold")


def calculate(data=DATA_PATH):
    data = Path(data)
    units = pd.read_parquet(data / "isomer_units.parquet")
    repeats = pd.read_parquet(data / "repeat_spectra.parquet")
    spectra = pd.read_parquet(data / "spectra.parquet")
    molecules = pd.read_parquet(data / "molecules.parquet")
    validate_tables(units, repeats, spectra, molecules)
    calibration = repeats.loc[repeats.fold == "train", "score"].to_numpy()
    threshold = float(np.quantile(calibration, QUANTILE, method="linear"))
    units["within"] = units.score >= threshold
    repeats["within"] = repeats.score >= threshold
    units["pair_first"] = np.minimum(units.query_molecule.to_numpy(), units.reference_molecule.to_numpy())
    units["pair_second"] = np.maximum(units.query_molecule.to_numpy(), units.reference_molecule.to_numpy())
    units["fold_pair"] = ["-".join(sorted(pair)) for pair in zip(units.query_fold, units.reference_fold)]
    known = spectra.instrument_type.notna() & spectra.collision_energy.notna()
    duplicates = units.identical_peak_list
    within = units.within
    no_duplicate_within = int((within & ~duplicates).sum())
    result = {
        "input": {"dataset": "MassSpecGym 1.5", "url": URL, "revision": REVISION,
                  "sha256": SHA256, "bytes": INPUT_BYTES},
        "method": {"repeat_quantile": QUANTILE, "intensity_power": INTENSITY_POWER,
                   "tolerance_da": TOLERANCE_DA, "tolerance_ppm": TOLERANCE_PPM,
                   "precursor_margin_da": PRECURSOR_MARGIN_DA,
                   "duplicate_atol": DUPLICATE_ATOL, "duplicate_rtol": DUPLICATE_RTOL},
        "coverage": {"recorded_condition": share(known.sum(), len(spectra)),
                     "adducts": sorted(spectra.adduct.dropna().unique().tolist()),
                     "instruments": sorted(spectra.instrument_type.dropna().unique().tolist()),
                     "molecule_keys": int(molecules.molecule_key.nunique()),
                     "molecule_records": len(molecules),
                     "keys_with_multiple_formulas": int((molecules.groupby("molecule_key").formula.nunique() > 1).sum()),
                     "keys_in_multiple_folds": int((molecules.groupby("molecule_key").fold.nunique() > 1).sum())},
        "threshold": threshold,
        "repeats": grouped_shares(repeats, "fold"),
        "units": share(within.sum(), len(units)),
        "query_spectra": int(units.query_id.nunique()),
        "molecules": int(pd.concat([units.query_molecule, units.reference_molecule]).nunique()),
        "formulas": int(units.formula.nunique()),
        "unscorable_repeats": int((repeats.score < 0).sum()),
        "unscorable_units": int((units.score < 0).sum()),
        "pairs": pair_summary(units, ["pair_first", "pair_second"]),
        "ordered_pairs": pair_summary(units, ["query_molecule", "reference_molecule"]),
        "breakdowns": {
            "instrument": grouped_shares(units, "instrument_type"),
            "adduct": grouped_shares(units, "adduct"),
            "query_fold": grouped_shares(units, "query_fold"),
            "fold_pair": grouped_shares(units, "fold_pair"),
            "same_fold": share(units.loc[units.query_fold == units.reference_fold, "within"].sum(),
                               (units.query_fold == units.reference_fold).sum()),
            "cross_fold": share(units.loc[units.query_fold != units.reference_fold, "within"].sum(),
                                (units.query_fold != units.reference_fold).sum()),
        },
        "duplicates": {
            "of_all_units": share(duplicates.sum(), len(units)),
            "of_within_units": share((duplicates & within).sum(), within.sum()),
            "molecule_pairs": len(units.loc[duplicates, ["pair_first", "pair_second"]].drop_duplicates()),
            "exact_normalised_peak_lists": int(units.exact_peak_list.sum()),
            "tolerance_only_peak_lists": int((duplicates & ~units.exact_peak_list).sum()),
            "without_duplicates_all_denominator": share(no_duplicate_within, len(units)),
            "without_duplicates_reduced_denominator": share(no_duplicate_within, (~duplicates).sum()),
            "nonidentical_score_at_least_0_99": int(((units.score >= 0.99) & ~duplicates).sum()),
        },
        "best_reference_different_condition": share((within & ~units.reference_same_condition).sum(), within.sum()),
        "query_repeat_agreement": agreement(units),
    }
    bands = [
        ("collision_energy", "collision_energy", [0, 20, 30, 40, 60, np.inf],
         ["(0,20]", "(20,30]", "(30,40]", "(40,60]", "over 60"]),
        ("fragment_peaks", "fragment_peak_count", [-1, 2, 5, 10, 20, 50, np.inf],
         ["0 to 2", "3 to 5", "6 to 10", "11 to 20", "21 to 50", "over 50"]),
        ("precursor_mz", "precursor_mz", [0, 200, 300, 400, 500, np.inf],
         ["up to 200", "(200,300]", "(300,400]", "(400,500]", "over 500"]),
    ]
    for name, column, edges, labels in bands:
        categories = pd.cut(units[column], edges, labels=labels)
        result["breakdowns"][name] = grouped_shares(units, categories)
        result["breakdowns"][f"{name}_outside_bands"] = share(
            units.loc[categories.isna(), "within"].sum(), categories.isna().sum()
        )
    result["collision_energy_zero_units"] = int((units.collision_energy == 0).sum())
    result["sensitivity"] = {}
    for quantile in [0.005, 0.01, 0.02, 0.05]:
        cutoff = float(np.quantile(calibration, quantile, method="linear"))
        result["sensitivity"][f"{100 * quantile:g}%"] = {
            "threshold": cutoff, **share((units.score >= cutoff).sum(), len(units))
        }
    example = units[(units.query_id == EXAMPLE_IDS[0]) & (units.reference_id == EXAMPLE_IDS[1])]
    if len(example) != 1:
        raise ValueError("The named example must be a unique query and best-reference pair")
    example_meta = pd.read_csv(data / "example_spectra.csv", float_precision="round_trip")
    example_peaks = pd.read_csv(data / "example_peaks.csv", float_precision="round_trip")
    example_pairs = pd.read_csv(data / "example_pairs.csv", float_precision="round_trip")
    validate_examples(example_pairs, example_meta, example_peaks, units, spectra)
    result["examples"] = example_pairs.set_index("pair_id").to_dict(orient="index")
    result["example"] = {
        "query_id": EXAMPLE_IDS[0], "reference_id": EXAMPLE_IDS[1],
        "score": float(example.score.iloc[0]), "formula": str(example.formula.iloc[0]),
        "condition": str(example.condition.iloc[0]),
        "query_fragment_peaks": int(example.fragment_peak_count.iloc[0]),
        "identical_peak_list": bool(example.identical_peak_list.iloc[0]),
    }
    return result


def check_expected(result, path):
    """Compare with fixed reference values, respecting their stated precision."""
    expected = json.loads(Path(path).read_text())
    checks = []
    for entry in expected["checks"]:
        actual = result
        for key in entry["path"]:
            actual = actual[key]
        value = entry["value"]
        if "decimals" in entry:
            matches = round(actual, entry["decimals"]) == value
        elif "absolute_tolerance" in entry:
            matches = abs(actual - value) <= entry["absolute_tolerance"]
        else:
            matches = actual == value
        checks.append({**entry, "actual": actual, "matches": bool(matches)})
    return {"matches_expected": all(item["matches"] for item in checks), "checks": checks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DATA_PATH, help="Included or recomputed tables")
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "summary.json", help="Results JSON")
    parser.add_argument("--expected", type=Path, default=DATA_PATH / "expected.json", help="Fixed reference values")
    args = parser.parse_args()
    result = calculate(args.data)
    result["validation"] = check_expected(result, args.expected)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"{result['units']['within']:,} / {result['units']['n']:,} units: {result['units']['percent']:.3f}%")
    print(f"Threshold: {result['threshold']:.16f}")
    passed = result["validation"]["matches_expected"]
    print(f"Expected values: {'MATCH' if passed else 'MISMATCH'}; results: {args.output}")
    if not passed:
        for check in result["validation"]["checks"]:
            if not check["matches"]:
                print(f"  {'/'.join(check['path'])}: expected {check['value']}, found {check['actual']}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
