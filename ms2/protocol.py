"""Fixed input and measurement choices used by the published tables."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REVISION = "018a9542378725cf5620423f930608730c36dd8c"
URL = (
    "https://huggingface.co/datasets/roman-bushuiev/MassSpecGym/resolve/"
    f"{REVISION}/data/MassSpecGym1.5.tsv"
)
SHA256 = "50cfdd1d22f79543c59555f9ce43c6893bd788a19a42b21fb4e2e3a54673c06a"
INPUT_BYTES = 261_197_365
INPUT_PATH = ROOT / "input" / "MassSpecGym1.5.tsv"
DATA_PATH = ROOT / "data"
QUANTILE = 0.01
INTENSITY_POWER = 0.5
TOLERANCE_DA = 0.01
TOLERANCE_PPM = 10.0
PRECURSOR_MARGIN_DA = 0.5
DUPLICATE_ATOL = 1e-6
DUPLICATE_RTOL = 1e-5
EXAMPLE_IDS = ("MassSpecGymID0052678", "MassSpecGymID0045860")
EXAMPLE_NAMES = ("Propazine", "Terbuthylazine")
EXAMPLES = (
    ("triazines", *EXAMPLE_IDS, *EXAMPLE_NAMES),
    ("sulfones", "MassSpecGymID0030777", "MassSpecGymID0081201",
     "4,4'-Diaminodiphenyl sulfone", "3,3'-Diaminodiphenyl sulfone"),
    ("flavonoids", "MassSpecGymID0171929", "MassSpecGymID0050074",
     "Naringenin chalcone", "Naringenin"),
)
