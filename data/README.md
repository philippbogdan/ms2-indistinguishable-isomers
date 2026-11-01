# Find and inspect the ambiguous comparisons

Start with [isomer_units.parquet](isomer_units.parquet): it identifies every one of the 62,535 comparisons, its best reference and its score. Apply the cutoff regenerated from [repeat_spectra.parquet](repeat_spectra.parquet) to recover the 11,554 flagged units. Join spectrum identifiers to the metadata for structures, conditions and folds. Complete peak lists for all three visual examples are bundled; the rest are available through the public input.

| File | Contents and join keys |
|:--|:--|
| `isomer_units.parquet` | One row per `(query_id, reference_molecule)`. `reference_id` is the first best-scoring reference in input order; `same_condition_reference_id` witnesses eligibility. |
| `repeat_spectra.parquet` | One row per repeat `query_id`, with the best other `reference_id` and a distinct `same_condition_repeat_id`. |
| `spectra.parquet` | Metadata for all 231,104 input spectra. Unique `spectrum_id`; `source_row` is the zero-based TSV data-row index. Includes keys, SMILES, formula, condition, fold, masses and peak counts. |
| `molecules.parquet` | Distinct `(molecule_key, smiles, formula, fold)` records. Keys are connectivity blocks. A key can have multiple SMILES records, so do not assume one row per key. Formula records are preserved as supplied. |
| `example_pairs.csv` | Three selected pairs in display order, with query/reference identifiers, names, scores, conditions, folds, retained peak counts and duplicate flags. `pair_id` joins to the example metadata. |
| `example_spectra.csv` | Six spectra with `pair_id`, `role`, `name`, public metadata and SMILES. `role` gives the query or reference position in the mirror plot. |
| `example_peaks.csv` | `spectrum_id`, zero-based `peak_index`, `mz`, `relative_intensity`, `retained`. Complete sorted lists; `retained` marks peaks used for scoring. |
| `expected.json` | Fixed reference counts and rounded values against which the offline command checks its independently calculated results. |

The unit table includes `query_molecule`, `reference_molecule`, `formula`, `condition`, `adduct`, `instrument_type`, `collision_energy`, `precursor_mz`, `fragment_peak_count`, `query_fold`, `reference_fold`, `score`, `reference_same_condition` and the duplicate flags. Instrument, energy, mass and fragment count refer to the query. Join `reference_id` to `spectra.spectrum_id` for the best reference's metadata. The repeat table records its `molecule_key`, query `condition`, `adduct`, `fold`, `score` and `reference_same_condition`.

The examples illustrate triazines, diaminodiphenyl sulfones and naringenin/chalcone. They were selected for their different structures and nearly matching spectra, not to estimate prevalence. Every pair shares its formula and recorded condition, has at least eight retained peaks on each side, and has unequal full peak lists under the duplicate tolerance. `ms2.results` checks their metadata and peak masks against the evidence; the small tests rescore them from their peaks.

`condition` serialises adduct, instrument type and energy with `|` separators. In the complete metadata table only, missing instruments become `NA` and missing energies become `-1.0` in this label; the corresponding typed columns retain missing values. Such rows can be references but cannot be eligible queries. Similarities are float64 and are not rounded in storage.

Recomputation can differ in the last floating-point bits between platforms. An effectively tied best reference can also change. The [cross-platform check](../REPRODUCE.md#verification) found no change to any reported number, count or verdict.

`identical_peak_list` uses equal-length, sorted, max-normalised full lists with numerical tolerance: absolute 1e-6 plus relative 1e-5 for both masses and intensities. It flags 897 units. `exact_peak_list` uses exact array equality and flags 889. These are separate from the similarity, which removes the precursor region. [REPRODUCE.md](../REPRODUCE.md) defines the arithmetic and ties.

For example, from the repository root:

```python
import json
import pandas as pd

units = pd.read_parquet("data/isomer_units.parquet")
threshold = json.load(open("results/summary.json"))["threshold"]
overlap = units.loc[units.score >= threshold]
spectra = pd.read_parquet("data/spectra.parquet").set_index("spectrum_id")
example = spectra.loc["MassSpecGymID0052678"]
```

Derived from R. Bushuiev and colleagues' [MassSpecGym](https://huggingface.co/datasets/roman-bushuiev/MassSpecGym), revision `018a9542378725cf5620423f930608730c36dd8c`. Keep the adjacent [MIT notice](MassSpecGym-LICENSE.txt) when redistributing these tables. Our derived contributions are dedicated under [CC0 1.0](../LICENSE).
