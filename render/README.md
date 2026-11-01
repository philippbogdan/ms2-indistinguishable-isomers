# Regenerate the four figures

Run the [figure commands](../REPRODUCE.md#regenerate-the-figures) from the repository root. Rendering needs the packages in [requirements-figures.txt](../requirements-figures.txt), including RDKit for the molecular structures; these are separate from the measurement dependencies.

| Script | Input and output |
|:--|:--|
| [make_data.py](make_data.py) | Validates `data/` against `data/expected.json`, then exports eight CSVs to `figures/data/`. |
| [fig1_atlas.py](fig1_atlas.py) | All isomer and repeat scores, plus the cutoff, to `figures/isomer-atlas.png`. |
| [fig2_examples.py](fig2_examples.py) | Example identities, SMILES and peaks, checked against isomer scores, to `figures/isomer-examples.png` and `.svg`. |
| [fig3_breakdown.py](fig3_breakdown.py) | Category counts independently recounted from isomer scores, to `figures/overlap-by-measurement.png` and `.svg`. |
| [fig4_sensitivity.py](fig4_sensitivity.py) | All isomer scores and four repeat quantiles to `figures/threshold-sensitivity.png`. |
| [common.py](common.py) | Shared style, empirical step curves, export and numerical/layout audit helpers. |

The eight generated tables are `isomer_scores.csv`, `repeat_scores.csv`, `threshold.csv`, `example_pairs.csv`, `example_spectra.csv`, `example_peaks.csv`, `breakdowns.csv` and `sensitivity.csv`. They are disposable exports of the published evidence, not a second source of results. Floating-point values are read back without loss of precision.

The five commands accept `--help`. `make_data.py` accepts `--data`, `--expected` and `--output` for checking another set of tables; renderers always read `figures/data/`. The example and breakdown renderers accept `--formats png`, `--formats svg` or both. The atlas and sensitivity renderers produce PNG only.

Rendering replaces the selected images and writes checks to `figures/audit/`. The atlas checks total and flagged counts and the repeat curves; the examples check identities, peaks and unchanged molecular connectivity; the breakdown checks every category; the sensitivity plot checks each marked cutoff against the empirical curve. All renderers check permitted text and canvas bounds. [Figure details](../figures/README.md) document scaling, denominators and display limits.
