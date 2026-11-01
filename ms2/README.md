# Apply the same ambiguity test to another library

Reuse the fragment matcher and repeat calibration to find isomer comparisons that resemble repeated measurements. The code separates scoring from enumeration and reporting, so each choice can be inspected or replaced. [Exact definitions and commands](../REPRODUCE.md).

| Module | What it provides |
|:--|:--|
| `similarity.py` | Precursor removal, greedy one-to-one matching, square-root cosine and maximum-reference selection. |
| `recompute.py` | Eligible repeat queries and directed isomer comparisons, scored on CPU from public input. |
| `tables.py` | Joinable scores, condition witnesses, molecule folds and duplicate flags. |
| `examples.py` | Three fixed structure pairs, complete peak extraction and consistency checks for the visual evidence. |
| `results.py` | Offline regeneration of reported counts, shares, cutoffs and example scores, checked against fixed expected values. |
| `protocol.py` | Pinned input, hash, measurement choices and example identities. |
| `download.py` | Verified download, moved into place only after its length and hash match. |

When porting, preserve query direction: the query peak mass determines matching tolerance. Reference maxima take the first attaining spectrum in input order. Eligibility requires a same-condition witness, but the best-reference search uses every spectrum of the reference molecule with the same adduct. Repeat calibration must use that same search rule. The [tests](../tests/README.md) exercise these consequential choices.
