"""Recompute the included tables from the hash-verified MassSpecGym TSV."""

import argparse
from itertools import permutations
from pathlib import Path
import time

import numpy as np
import pandas as pd

from .download import verify
from .examples import extract_examples
from .protocol import DATA_PATH, INPUT_PATH, PRECURSOR_MARGIN_DA
from .similarity import maximum_scores
from .tables import write_tables


def load_input(path):
    """Read peaks in double precision, sort stably and normalise each spectrum."""
    frame = pd.read_csv(path, sep="\t")
    mass_parts = frame.mzs.str.split(",")
    intensity_parts = frame.intensities.str.split(",")
    counts = mass_parts.str.len().to_numpy()
    if not np.array_equal(counts, intensity_parts.str.len().to_numpy()):
        raise ValueError("Mass and intensity lists have different lengths")
    offsets = np.zeros(len(frame) + 1, dtype=np.int64)
    offsets[1:] = np.cumsum(counts)
    mz = np.fromiter((float(v) for part in mass_parts for v in part), np.float64, offsets[-1])
    intensity = np.fromiter(
        (float(v) for part in intensity_parts for v in part), np.float64, offsets[-1]
    )
    for row in range(len(frame)):
        start, stop = offsets[row], offsets[row + 1]
        order = np.argsort(mz[start:stop], kind="stable")
        mz[start:stop] = mz[start:stop][order]
        values = intensity[start:stop][order]
        if len(values) == 0 or values.max() <= 0:
            raise ValueError(f"Spectrum row {row} has no positive intensity")
        intensity[start:stop] = values / values.max()
    meta = frame.drop(columns=["mzs", "intensities"]).copy()
    meta["row"] = np.arange(len(meta))
    meta["npeaks"] = counts
    meta["cond"] = (
        meta.adduct + "|" + meta.instrument_type.fillna("NA") + "|"
        + meta.collision_energy.fillna(-1.0).map(repr)
    )
    precursor = meta.precursor_mz.to_numpy()
    meta["nfrag"] = [
        int((mz[offsets[i]:offsets[i + 1]] <= precursor[i] - PRECURSOR_MARGIN_DA).sum())
        for i in range(len(meta))
    ]
    return meta, mz, intensity, offsets


def comparison_sets(meta):
    """Define eligible queries before broadening reference searches by adduct."""
    references = meta.groupby(["inchikey", "adduct"]).row.apply(np.array).to_dict()
    known = meta[meta.instrument_type.notna() & meta.collision_energy.notna()]
    repeat_sizes = known.groupby(["inchikey", "cond"]).row.transform("size")
    repeats = known.loc[repeat_sizes >= 2, "row"].to_numpy()
    rows = []
    for (formula, condition), group in known.groupby(["formula", "cond"]):
        molecules = sorted(group.inchikey.unique())
        for first, second in permutations(molecules, 2):
            for query in group.loc[group.inchikey == first, "row"]:
                rows.append((query, first, second, condition, formula, group.adduct.iloc[0]))
    units = pd.DataFrame(rows, columns=["row", "A", "B", "cond", "formula", "adduct"])
    return references, repeats, units


def recompute(path, output):
    started = time.perf_counter()
    verify(path)
    meta, mz, intensity, offsets = load_input(path)
    print(f"Read {len(meta):,} spectra", flush=True)
    references, repeat_rows, units = comparison_sets(meta)
    keys, adducts = meta.inchikey.to_numpy(), meta.adduct.to_numpy()
    precursor = meta.precursor_mz.to_numpy()
    repeat_sets = []
    for row in repeat_rows:
        candidates = references[(keys[row], adducts[row])]
        repeat_sets.append(candidates[candidates != row])
    repeat_score, repeat_reference = maximum_scores(
        mz, intensity, offsets, precursor, repeat_rows, repeat_sets
    )
    repeats = pd.DataFrame({
        "row": repeat_rows, "s_repeat": repeat_score, "best_ref": repeat_reference
    })
    units["s_meas"], units["best_ref"] = maximum_scores(
        mz, intensity, offsets, precursor, units.row.to_numpy(),
        [references[(key, adduct)] for key, adduct in zip(units.B, units.adduct)],
    )
    if (repeat_score < 0).any() or (units.s_meas < 0).any():
        raise ValueError("An unscorable spectrum was found; the published data have none")
    write_tables(meta, mz, intensity, offsets, units, repeats, output)
    print(f"Wrote tables to {output} in {time.perf_counter() - started:.1f} s")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT_PATH, help="Pinned TSV input")
    parser.add_argument("--output", type=Path, default=DATA_PATH, help="Derived table directory")
    parser.add_argument("--examples-only", action="store_true",
                        help="Extract the selected peaks using existing tables; do not rescore the input")
    args = parser.parse_args()
    if args.examples_only:
        verify(args.input)
        extract_examples(args.input, args.output)
        print(f"Wrote selected example tables to {args.output}")
    else:
        recompute(args.input, args.output)


if __name__ == "__main__":
    main()
