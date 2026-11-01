"""Preserve three selected structure pairs and their complete measured peak lists."""

from pathlib import Path

import numpy as np
import pandas as pd

from .protocol import DUPLICATE_ATOL, DUPLICATE_RTOL, EXAMPLES, PRECURSOR_MARGIN_DA


def select_examples(units, spectra):
    """Select fixed illustrations, not a sample for estimating prevalence."""
    metadata = spectra.set_index("spectrum_id", drop=False)
    pair_rows, spectrum_rows = [], []
    for pair_id, query_id, reference_id, query_name, reference_name in EXAMPLES:
        matches = units[(units.query_id == query_id) & (units.reference_id == reference_id)]
        if len(matches) != 1:
            raise ValueError(f"Example {pair_id} must identify one scored unit")
        unit = matches.iloc[0]
        query, reference = metadata.loc[query_id], metadata.loc[reference_id]
        if (unit.identical_peak_list or not unit.reference_same_condition
                or query.condition != reference.condition or query.formula != reference.formula
                or query.molecule_key == reference.molecule_key
                or min(query.fragment_peak_count, reference.fragment_peak_count) < 8):
            raise ValueError(f"Example {pair_id} does not meet the illustration rule")
        pair_rows.append({
            "pair_id": pair_id, "query_id": query_id, "reference_id": reference_id,
            "query_name": query_name, "reference_name": reference_name,
            "formula": unit.formula, "condition": unit.condition,
            "adduct": query.adduct, "instrument_type": query.instrument_type,
            "collision_energy": float(query.collision_energy),
            "query_fold": query.fold, "reference_fold": reference.fold,
            "score": float(unit.score),
            "query_fragment_peaks": int(query.fragment_peak_count),
            "reference_fragment_peaks": int(reference.fragment_peak_count),
            "identical_peak_list": bool(unit.identical_peak_list),
            "exact_peak_list": bool(unit.exact_peak_list),
        })
        for role, name, row in [("query", query_name, query), ("reference", reference_name, reference)]:
            spectrum_rows.append({"pair_id": pair_id, "role": role, "name": name, **row.to_dict()})
    return pd.DataFrame(pair_rows), pd.DataFrame(spectrum_rows)


def validate_examples(pairs, example_spectra, peaks, units, spectra):
    """Check identities, counts, normalisation and duplicate flags before display."""
    expected_pairs, expected_spectra = select_examples(units, spectra)
    for actual, expected in [(pairs, expected_pairs), (example_spectra, expected_spectra)]:
        pd.testing.assert_frame_equal(actual.reset_index(drop=True), expected,
                                      check_dtype=False, check_exact=True)
    if set(peaks.spectrum_id) != set(example_spectra.spectrum_id):
        raise ValueError("Example peaks and metadata contain different spectra")
    lists = {}
    for row in example_spectra.itertuples(index=False):
        part = peaks[peaks.spectrum_id == row.spectrum_id].sort_values("peak_index")
        if (len(part) != row.peak_count or int(part.retained.sum()) != row.fragment_peak_count
                or not np.array_equal(part.peak_index, np.arange(len(part)))
                or not np.array_equal(part.retained, part.mz <= row.precursor_mz - PRECURSOR_MARGIN_DA)
                or not np.isfinite(part[["mz", "relative_intensity"]]).all().all()
                or (part.mz.diff().dropna() < 0).any()
                or (part.relative_intensity < 0).any() or part.relative_intensity.max() != 1):
            raise ValueError(f"Invalid example peak list: {row.spectrum_id}")
        lists[row.spectrum_id] = part[["mz", "relative_intensity"]].to_numpy()
    for row in pairs.itertuples(index=False):
        first, second = lists[row.query_id], lists[row.reference_id]
        if first.shape == second.shape and np.allclose(
                first, second, atol=DUPLICATE_ATOL, rtol=DUPLICATE_RTOL):
            raise ValueError(f"Example {row.pair_id} has numerically identical peak lists")


def write_examples(units, spectra, peak_lists, output):
    """Write only the fixed examples, preserving full lists and the scoring mask."""
    pairs, example_spectra = select_examples(units, spectra)
    rows = []
    for row in example_spectra.itertuples(index=False):
        mz, intensity = peak_lists[row.spectrum_id]
        rows.extend((row.spectrum_id, i, float(mass), float(value),
                     bool(mass <= row.precursor_mz - PRECURSOR_MARGIN_DA))
                    for i, (mass, value) in enumerate(zip(mz, intensity)))
    peaks = pd.DataFrame(rows, columns=[
        "spectrum_id", "peak_index", "mz", "relative_intensity", "retained",
    ])
    validate_examples(pairs, example_spectra, peaks, units, spectra)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(output / "example_pairs.csv", index=False)
    example_spectra.to_csv(output / "example_spectra.csv", index=False)
    peaks.to_csv(output / "example_peaks.csv", index=False)


def extract_examples(path, output):
    """Extract selected peaks from a verified input using the existing score tables."""
    output = Path(output)
    units = pd.read_parquet(output / "isomer_units.parquet")
    spectra = pd.read_parquet(output / "spectra.parquet")
    _, selected = select_examples(units, spectra)
    wanted, peak_lists = set(selected.spectrum_id), {}
    with pd.read_csv(path, sep="\t", usecols=["identifier", "mzs", "intensities"],
                     chunksize=10_000) as chunks:
        for chunk in chunks:
            for row in chunk[chunk.identifier.isin(wanted)].itertuples(index=False):
                mz = np.fromstring(row.mzs, sep=",", dtype=np.float64)
                intensity = np.fromstring(row.intensities, sep=",", dtype=np.float64)
                if len(mz) != len(intensity) or not len(mz) or intensity.max() <= 0:
                    raise ValueError(f"Invalid input peaks: {row.identifier}")
                order = np.argsort(mz, kind="stable")
                peak_lists[row.identifier] = (mz[order], intensity[order] / intensity.max())
            if wanted <= peak_lists.keys():
                break
    if wanted != peak_lists.keys():
        raise ValueError("The pinned input is missing an example spectrum")
    write_examples(units, spectra, peak_lists, output)
