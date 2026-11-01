"""Write compact, joinable evidence tables without changing scored values."""

from pathlib import Path

import numpy as np
import pandas as pd

from .examples import select_examples, write_examples
from .protocol import DUPLICATE_ATOL, DUPLICATE_RTOL


def peak_list_flags(mz, intensity, offsets, queries, references):
    """Compare complete sorted, max-normalised lists, before precursor removal."""
    similar = np.zeros(len(queries), dtype=bool)
    exact = np.zeros(len(queries), dtype=bool)
    for k, (query, reference) in enumerate(zip(queries, references)):
        a, b = offsets[query], offsets[query + 1]
        c, d = offsets[reference], offsets[reference + 1]
        if b - a != d - c:
            continue
        exact[k] = np.array_equal(mz[a:b], mz[c:d]) and np.array_equal(intensity[a:b], intensity[c:d])
        similar[k] = exact[k] or (
            np.allclose(mz[a:b], mz[c:d], atol=DUPLICATE_ATOL, rtol=DUPLICATE_RTOL)
            and np.allclose(intensity[a:b], intensity[c:d], atol=DUPLICATE_ATOL, rtol=DUPLICATE_RTOL)
        )
    return similar, exact


def write_tables(meta, mz, intensity, offsets, units, repeats, output):
    """Materialise scores, metadata and the example spectra in a public schema."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    fields = {
        "row": "source_row", "identifier": "spectrum_id", "inchikey": "molecule_key",
        "smiles": "smiles", "formula": "formula", "adduct": "adduct",
        "instrument_type": "instrument_type", "collision_energy": "collision_energy",
        "cond": "condition", "fold": "fold", "precursor_mz": "precursor_mz",
        "parent_mass": "parent_mass", "npeaks": "peak_count", "nfrag": "fragment_peak_count",
    }
    spectra = meta[list(fields)].rename(columns=fields)
    if not np.array_equal(spectra.source_row, np.arange(len(spectra))):
        raise ValueError("Spectrum rows must retain their original input order")
    if not spectra.spectrum_id.is_unique:
        raise ValueError("Spectrum identifiers must be unique")
    spectra.to_parquet(output / "spectra.parquet", index=False)
    molecules = spectra[["molecule_key", "smiles", "formula", "fold"]].drop_duplicates()
    molecules.sort_values(["molecule_key", "formula", "smiles"]).to_parquet(
        output / "molecules.parquet", index=False
    )
    known = meta[meta.instrument_type.notna() & meta.collision_energy.notna()]
    witnesses = known.groupby(["inchikey", "formula", "cond"]).row.first().to_dict()
    repeat_groups = known.groupby(["inchikey", "cond"]).row.apply(np.array).to_dict()
    identifiers = meta.identifier.to_numpy()
    unit_queries = units.row.to_numpy(dtype=np.int64)
    unit_references = units.best_ref.to_numpy(dtype=np.int64)
    query = spectra.iloc[unit_queries].reset_index(drop=True)
    reference = spectra.iloc[unit_references].reset_index(drop=True)
    public_units = pd.DataFrame({
        "query_id": query.spectrum_id, "reference_id": reference.spectrum_id,
        "query_molecule": units.A.to_numpy(), "reference_molecule": units.B.to_numpy(),
        "formula": units.formula.to_numpy(), "condition": units.cond.to_numpy(),
        "adduct": query.adduct, "instrument_type": query.instrument_type,
        "collision_energy": query.collision_energy, "precursor_mz": query.precursor_mz,
        "fragment_peak_count": query.fragment_peak_count, "query_fold": query.fold,
        "reference_fold": reference.fold, "score": units.s_meas.to_numpy(),
        "reference_same_condition": query.condition == reference.condition,
        "same_condition_reference_id": [
            identifiers[witnesses[(key, formula, condition)]]
            for key, formula, condition in zip(units.B, units.formula, units.cond)
        ],
    })
    public_units["identical_peak_list"], public_units["exact_peak_list"] = peak_list_flags(
        mz, intensity, offsets, unit_queries, unit_references
    )
    public_units.to_parquet(output / "isomer_units.parquet", index=False)
    repeat_queries = repeats.row.to_numpy(dtype=np.int64)
    repeat_references = repeats.best_ref.to_numpy(dtype=np.int64)
    query = spectra.iloc[repeat_queries].reset_index(drop=True)
    reference = spectra.iloc[repeat_references].reset_index(drop=True)
    repeat_witnesses = []
    for row, key, condition in zip(repeat_queries, query.molecule_key, query.condition):
        candidates = repeat_groups[(key, condition)]
        repeat_witnesses.append(identifiers[candidates[1] if candidates[0] == row else candidates[0]])
    pd.DataFrame({
        "query_id": query.spectrum_id, "reference_id": reference.spectrum_id,
        "molecule_key": query.molecule_key, "condition": query.condition,
        "adduct": query.adduct, "fold": query.fold, "score": repeats.s_repeat.to_numpy(),
        "reference_same_condition": query.condition == reference.condition,
        "same_condition_repeat_id": repeat_witnesses,
    }).to_parquet(output / "repeat_spectra.parquet", index=False)
    _, example = select_examples(public_units, spectra)
    peak_lists = {}
    for row in example.itertuples(index=False):
        start, stop = offsets[row.source_row], offsets[row.source_row + 1]
        peak_lists[row.spectrum_id] = (mz[start:stop], intensity[start:stop])
    write_examples(public_units, spectra, peak_lists, output)
