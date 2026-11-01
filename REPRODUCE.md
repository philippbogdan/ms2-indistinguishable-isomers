# Reproduce the measurement

The included tables are sufficient to regenerate the results offline. Recomputing the spectrum scores additionally needs the pinned public TSV. No model, GPU or private input is needed.

```text
Pinned public TSV + SHA-256
            |
    sort and normalise peaks
            |
    repeats + isomer units
            |
     greedy peak matching
            |
       data/*.parquet
            |
     results/summary.json
```

## Install

Run all commands from the repository root. Python 3.14 was used for the local checks; the dependencies are pinned to the versions used to prepare the included tables.

```sh
python -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

NumPy, pandas and Numba use BSD-style licences; PyArrow uses Apache 2.0. They are installed dependencies, not bundled code. The installed NumPy distribution also contains components under its documented third-party licences. Plotting dependencies are separate and are only needed for [figure regeneration](#regenerate-the-figures).

## Check the included evidence

```sh
.venv/bin/python -m ms2.results
```

This reads only `data/`, checks joins and conditions, recomputes the cutoff, counts and shares, validates all three example pairs and their peak masks, and writes [results/summary.json](results/summary.json). It prints `MATCH` or `MISMATCH` against [data/expected.json](data/expected.json) and exits unsuccessfully on a mismatch. Integer counts must match exactly. Rounded percentages and example scores must match at their stated precision; the full cutoff has an absolute tolerance of 5e-13. The fixed expected values are never replaced by this command. It takes about a second in the preparation environment; it does not rescore the original peak lists.

## Fetch the exact public input

```sh
.venv/bin/python -m ms2.download
```

The script downloads [MassSpecGym 1.5](https://huggingface.co/datasets/roman-bushuiev/MassSpecGym/resolve/018a9542378725cf5620423f930608730c36dd8c/data/MassSpecGym1.5.tsv) to `input/MassSpecGym1.5.tsv`. The revision is `018a9542378725cf5620423f930608730c36dd8c`, the length is 261,197,365 bytes, and the SHA-256 is:

```text
50cfdd1d22f79543c59555f9ce43c6893bd788a19a42b21fb4e2e3a54673c06a
```

An existing destination is verified instead of downloaded again. A download is moved into place only after it passes both checks. A wrong existing file produces an error and is not overwritten.

## Recompute the scores and tables

```sh
.venv/bin/python -m ms2.recompute
.venv/bin/python -m ms2.results
```

The first command verifies the input hash again and replaces the derived Parquet and example CSV files in `data/`. It leaves the licence, documentation and fixed expected values intact. The second command regenerates and checks all reported results. Recomputation took **33.3 seconds** in the clean-Linux check described [below](#verification); runtime depends on the machine, and first use also compiles the Numba functions.

For a comparison without overwriting included tables:

```sh
.venv/bin/python -m ms2.recompute --output /tmp/ms2-recomputed
.venv/bin/python -m ms2.results --data /tmp/ms2-recomputed --output /tmp/ms2-summary.json
```

Each command also accepts `--help`. The download has `--output`; recomputation has `--input`, `--output` and `--examples-only`; results has `--data`, `--output` and `--expected`.

To extract only the selected examples from the verified input, keeping the existing score tables:

```sh
.venv/bin/python -m ms2.recompute --examples-only
.venv/bin/python -m ms2.results
```

This streams the input and replaces only `example_pairs.csv`, `example_spectra.csv` and `example_peaks.csv`. Full recomputation writes these same files. The selection in `protocol.py` fixes the three displayed pairs; `examples.py` checks that both members have at least eight retained peaks, the same formula and condition, distinct molecule keys and no tolerance-based duplicate flag. These illustrations are not a random sample or the highest-scoring pairs.

## Exact definitions

**Molecule and condition.** A molecule is the connectivity block in the input `inchikey` column. Stereochemistry does not create separate molecules. A condition is `(adduct, instrument_type, collision_energy)`. Only queries with recorded instrument type and energy are eligible. Energy is compared as recorded, without conversion. Reference searches can include spectra with missing condition metadata.

**Peak handling.** Parse comma-separated masses and intensities as float64. Sort peaks by mass with a stable sort and divide intensities by their maximum. For scoring, remove peaks strictly above that spectrum's own precursor m/z minus 0.5 Da. Retain a peak exactly at the cutoff. Renormalise retained intensities by their maximum, then take their square roots. Renormalisation cancels in the cosine but is retained to reproduce the arithmetic.

**Matching.** Enumerate candidate peak pairs in query-mass order, then reference-mass order. The inclusive tolerance is `max(0.01 Da, query_peak_mz * 10e-6)`. Sort candidates by absolute mass difference using NumPy/Numba's default `argsort`, then accept each if neither peak is already used. Do not maximise intensity agreement or shift peaks by precursor mass. The tolerance is defined from the first spectrum, so keep the query/reference direction. Equal-distance ties follow the implementation's default sort order; the dependency versions are pinned.

**Similarity.** Sum the products of the weights of matched peaks and divide by the product of the spectra's Euclidean weight norms. Unmatched peaks contribute to the denominator. A spectrum with no retained positive intensity gives no score. There are no such queries in the included repeat or unit tables. The port stops rather than silently assigning a verdict if one appears.

**Repeats and cutoff.** Include each spectrum whose molecule has another spectrum under exactly the same condition. Its repeat score is the maximum similarity to every other spectrum of that molecule with the same adduct, regardless of instrument or energy. Exclude the query itself. Set the cutoff to `numpy.quantile(training_repeat_scores, 0.01, method="linear")`. Evaluate the same cutoff on validation and test repeats.

**Isomer units.** For each formula and known condition containing distinct molecule keys, enumerate every ordered pair `(A, B)` and every query spectrum of A in that group. A unit's score is its maximum similarity to all spectra of B with the query's adduct. The best reference need not have the query's condition. If several references attain the same maximum, store the first in TSV row order. Count a unit as within the spread when `score >= cutoff`. The stored same-condition witness proves eligibility independently of the chosen best reference.

**Folds and weighting.** Use each molecule's actual fold, separately for query and reference. All folds contribute to the isomer table. The [MassSpecGym paper](https://arxiv.org/abs/2410.23326) describes a split with a minimum MCES edge edit distance of 10 across folds; this repository does not recompute that structural distance. Molecules with more query spectra receive greater unit weight. Pair summaries pool units across all eligible conditions and, for unordered pairs, both directions.

**Duplicates.** Compare the complete sorted, max-normalised peak lists, before precursor removal. The `identical_peak_list` flag reproduces `numpy.allclose` on both mass and intensity arrays with `atol=1e-6` and `rtol=1e-5`, requiring equal lengths. This flags 897 units; the name does not mean literal equality. `exact_peak_list` instead uses elementwise exact equality on both arrays and flags 889. Exclusion percentages use the original tolerance flag. Neither check establishes which library annotation is correct.

**Repeat agreement.** Within each `(A, B, condition)` group, count ordered pairs of distinct query spectra. If a group has `n` units and `k` within the spread, it contributes `k*(n-1)` conditioning cases and `k*(k-1)` agreeing cases for an initially within-spread query. Replace `k` by `n-k` for an initially outside-spread query. No independence assumption or confidence interval is attached to these descriptive shares.

**Breakdowns.** Energy and precursor-mass intervals are right-closed. Energy zero lies outside the positive-energy bands and is reported separately. Peak-count bins use the number of peaks left after precursor removal. Nearby quantiles use the same training-repeat scores and the same linear interpolation. The results retain both counts and fractions for every denominator.

## Small checks

```sh
.venv/bin/python -m unittest discover -s tests -v
```

The tests cover matching, precursor removal, reference selection, condition eligibility, duplicate tolerances, agreement counting, all three included examples, rejection of a corrupted example mask and download verification. They do not download the public file or run the full analysis. For a quick interpreted check without compilation, prefix the command with `NUMBA_DISABLE_JIT=1`; run normally as well when checking the compiled scoring path.

## Regenerate the figures

The four figures can be regenerated offline from the included tables. Install the additional plotting dependencies, export the checked drawing data, then run each renderer:

```sh
.venv/bin/python -m pip install -r requirements-figures.txt
.venv/bin/python render/make_data.py
.venv/bin/python render/fig1_atlas.py
.venv/bin/python render/fig2_examples.py
.venv/bin/python render/fig3_breakdown.py
.venv/bin/python render/fig4_sensitivity.py
```

`make_data.py` checks the evidence against `data/expected.json` before exporting eight CSV tables to `figures/data/`. The renderers read only those tables and replace the four PNGs and the example and breakdown SVGs in `figures/`. No public-input download is needed.

The scripts record numerical checks, drawing transformations and label bounds in `figures/audit/`. Generated CSVs, audit files and the Matplotlib cache are ignored by version control. [render/README.md](render/README.md) describes the scripts and options; [figures/README.md](figures/README.md) explains the encodings and display limits.

## Verification

The full path passed on a clean Linux machine in a fresh virtual environment installed from `requirements.txt`: the public-input hash check, recomputation in **33.3 seconds**, `ms2.results` on the recomputed tables (`MATCH`), `ms2.results` on the included tables (`MATCH`), and **12 tests** (`OK`). Separately, `ms2.download` fetched the public file over the network on macOS and verified its SHA-256.

The Linux tables match the included tables apart from floating-point details: **100 repeat scores** and **125 isomer-unit scores** differ by at most **3.3e-16**, and two recorded collision energies differ by one unit in the last place. One unit selects a different best reference from two references tied within **1e-16**. Every reported number, count and overlap verdict is identical. Reproduction therefore means agreement at the reported precision and on all verdicts, rather than byte-for-byte equality of floating-point tables or tied reference identifiers.
