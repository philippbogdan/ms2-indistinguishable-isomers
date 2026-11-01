"""Shared settings and helpers for the four figures.

Read round-trip floats from figures/data, render at 3600 pixels wide and
record checks in figures/audit. Text sizes are checked at 800 pixels wide.
"""

from __future__ import annotations

import json
import argparse
import os
from pathlib import Path
import tempfile

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"
PACK = FIGURES / "data"
AUDIT = FIGURES / "audit"
# Keep Matplotlib's cache within the repository on machines without a writable home.
os.environ.setdefault("MPLCONFIGDIR", str(FIGURES / ".matplotlib"))

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.text import Text  # noqa: E402

WIDTH_IN = 8.0
DPI = 450
DISPLAY_WIDTH_PX = 800
DISPLAY_SCALE = DISPLAY_WIDTH_PX / (WIDTH_IN * DPI)  # native px -> display px

# Each figure identifies colour through shape, line style or position as well.
BLUE = "#004488"  # Orbitrap (fig 1), query spectrum (fig 2), all comparisons (fig 4)
ORANGE = "#CC5500"  # QTOF (fig 1), reference spectrum (fig 2)
CHARCOAL = "#333333"  # numerically equal peak lists (figs 1 and 4)
GREEN = "#117733"  # repeat spectra (fig 1)
BLACK = "#000000"  # cutoff line (fig 1), category points (fig 3)
GREY = "#6E6E6E"  # overall share (fig 3), main cutoff rule (fig 4)
LIGHT_GREY = "#8C8C8C"  # other cutoff rules (fig 4)
GRID = "#E3E3E3"
PALETTE = {
    "blue": BLUE,
    "orange": ORANGE,
    "charcoal": CHARCOAL,
    "green": GREEN,
    "black": BLACK,
    "grey": GREY,
    "light_grey": LIGHT_GREY,
}

RC = {
    "font.family": "DejaVu Sans",
    "font.size": 12.5,
    "axes.labelsize": 13.5,
    "axes.titlesize": 13.5,
    "xtick.labelsize": 12.5,
    "ytick.labelsize": 12.5,
    "legend.fontsize": 12.5,
    "axes.linewidth": 0.8,
    "axes.edgecolor": "#222222",
    "axes.labelcolor": "#111111",
    "xtick.color": "#222222",
    "ytick.color": "#222222",
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 4.0,
    "ytick.major.size": 4.0,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "legend.frameon": False,
    "path.simplify": False,  # draw every empirical step vertex
    "svg.hashsalt": "ms2-indistinguishable-isomers",
    "axes.formatter.use_mathtext": False,
}

METADATA = {"Creator": f"Matplotlib {mpl.__version__}"}
RT = {"float_precision": "round_trip"}
SCORE_MAX = 1.0 + 1e-12


def read(name: str) -> pd.DataFrame:
    if not (PACK / name).is_file():
        raise FileNotFoundError("Figure data are missing; run python render/make_data.py first")
    return pd.read_csv(PACK / name, **RT)


def arguments(description: str, default_formats: tuple[str, ...]) -> list[str]:
    """Parse the shared output options before drawing a figure."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--formats", nargs="+", choices=default_formats,
                        default=list(default_formats), help="Image formats to replace in figures/")
    return parser.parse_args().formats


def check_scores(values: np.ndarray, label: str) -> None:
    values = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{label}: non-finite score")
    if values.min() < 0.0 or values.max() > SCORE_MAX:
        raise ValueError(f"{label}: score outside [0, 1 + 1e-12]")


def ge_counts(scores: np.ndarray, cutoffs: np.ndarray) -> np.ndarray:
    """count(score >= t) for each t, tie aware."""
    ordered = np.sort(np.asarray(scores, dtype=float))
    return ordered.size - np.searchsorted(ordered, cutoffs, side="left")


def survival_vertices(scores: np.ndarray, denominator: int, lo: float = 0.0, hi: float = 1.0):
    """Exact empirical step function t -> 100 * count(score >= t) / denominator.

    Steps sit at the union of the observed scores and the end points lo, hi;
    the function is left continuous, so at an observed score the count still
    includes it and the drop happens immediately above it. Returns stored
    cutoff coordinates (unclipped) and percentages; callers clip only the
    drawing coordinates.
    """
    ordered = np.sort(np.asarray(scores, dtype=float))
    unique = np.unique(ordered)
    unique = unique[unique >= lo]
    count_ge = ordered.size - np.searchsorted(ordered, unique, side="left")
    count_gt = ordered.size - np.searchsorted(ordered, unique, side="right")
    start = ordered.size - np.searchsorted(ordered, lo, side="left")
    t = np.concatenate([[lo], np.repeat(unique, 2), [max(hi, unique[-1] if unique.size else hi)]])
    c = np.concatenate([[start], np.column_stack([count_ge, count_gt]).ravel(), [count_gt[-1] if unique.size else start]])
    pct = 100.0 * c.astype(float) / denominator
    return t.astype(float), pct


def evaluate_survival(t_vertices: np.ndarray, pct_vertices: np.ndarray, t: float) -> float:
    """Value of the left-continuous step at t from its vertex list (for checks)."""
    # first vertex pair at or beyond t: the value on (previous, t] is the level
    # that ends at t, which is the first vertex with t_vertex >= t.
    index = int(np.searchsorted(t_vertices, t, side="left"))
    return float(pct_vertices[min(index, len(pct_vertices) - 1)])


def text_inventory(fig) -> list[dict]:
    """Every visible, non-empty text artist with its size at 800 px display."""
    renderer = fig.canvas.get_renderer()
    undrawn = set()  # tick labels outside the view interval are never drawn
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            drawn = {id(t) for t in axis._update_ticks()}
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                if id(tick) not in drawn:
                    undrawn.update({id(tick.label1), id(tick.label2)})
    items = []
    for artist in fig.findobj(Text):
        if not artist.get_visible() or id(artist) in undrawn:
            continue
        content = artist.get_text()
        if not content.strip():
            continue
        box = artist.get_window_extent(renderer)
        if box.width == 0 and box.height == 0:
            continue
        em_px = artist.get_fontsize() * DPI / 72.0
        items.append(
            {
                "text": content,
                "fontsize_pt": artist.get_fontsize(),
                "em_display_px": round(em_px * DISPLAY_SCALE, 2),
                "bbox_display_px": [round(v * DISPLAY_SCALE, 1) for v in box.bounds],
                "inside_canvas": bool(
                    box.x0 >= -0.5
                    and box.y0 >= -0.5
                    and box.x1 <= fig.bbox.width + 0.5
                    and box.y1 <= fig.bbox.height + 0.5
                ),
            }
        )
    return items


def audit_text(fig, allowed: set[str], numeric_ok: bool = True) -> dict:
    items = text_inventory(fig)
    bad = []
    for item in items:
        content = " ".join(item["text"].split())  # a wrapped label keeps its words
        if content in allowed:
            continue
        stripped = content.replace("−", "-").replace(",", "")
        if numeric_ok:
            try:
                float(stripped)
                continue
            except ValueError:
                pass
        bad.append(content)
    return {
        "texts": items,
        "disallowed": bad,
        "min_em_display_px": min(i["em_display_px"] for i in items),
        "all_inside_canvas": all(i["inside_canvas"] for i in items),
    }


def flatten_png(path: Path) -> dict:
    """Matplotlib writes RGBA PNGs; ours are opaque, so drop the alpha channel
    after checking it is 255 everywhere (a lossless change)."""
    from PIL import Image
    from PIL.PngImagePlugin import PngInfo

    with Image.open(path) as image:
        image.load()
        info = {k: v for k, v in image.info.items() if isinstance(v, str)}
        alpha = np.asarray(image.getchannel("A")) if image.mode == "RGBA" else None
        if alpha is not None and alpha.min() != 255:
            raise ValueError(f"{path.name}: not opaque")
        rgb = image.convert("RGB")
    text = PngInfo()
    for key, value in info.items():
        text.add_text(key, value)
    rgb.save(path, format="PNG", pnginfo=text, optimize=True, dpi=(DPI, DPI))
    return {"alpha_min_before": int(alpha.min()) if alpha is not None else None, "mode_after": "RGB"}


def save(fig, stem: str, formats: list[str]) -> dict:
    """Replace each image after rendering it fully, with opaque PNGs and SVG glyph paths."""
    FIGURES.mkdir(parents=True, exist_ok=True)
    outputs = []
    with mpl.rc_context({"svg.fonttype": "path", "savefig.bbox": None}):
        for kind in dict.fromkeys(formats):
            if kind not in {"png", "svg"}:
                raise ValueError(f"Unsupported figure format: {kind}")
            destination = FIGURES / f"{stem}.{kind}"
            fd, name = tempfile.mkstemp(prefix=f".{stem}-", suffix=f".{kind}", dir=FIGURES)
            os.close(fd)
            temporary = Path(name)
            try:
                metadata = dict(METADATA)
                if kind == "svg":
                    metadata["Date"] = None
                fig.savefig(temporary, format=kind, dpi=DPI, bbox_inches=None,
                            facecolor="white", edgecolor="none", transparent=False, metadata=metadata)
                if kind == "png":
                    flatten_png(temporary)
                temporary.replace(destination)
            finally:
                temporary.unlink(missing_ok=True)
            outputs.append({"path": f"figures/{destination.name}", "bytes": destination.stat().st_size})
    return {"matplotlib": mpl.__version__, "dpi": DPI, "outputs": outputs}


def write_audit(stem: str, audit: dict) -> None:
    """Store the numerical and layout checks alongside the generated figure data."""
    AUDIT.mkdir(parents=True, exist_ok=True)
    (AUDIT / f"{stem}.json").write_text(json.dumps(audit, indent=2, default=str) + "\n")


def export(fig, stem: str, formats: list[str], provenance: dict, audit: dict) -> dict:
    if audit["text"]["disallowed"] or not audit["text"]["all_inside_canvas"]:
        raise ValueError(f"{stem}: unexpected or clipped figure text")
    report = save(fig, stem, formats)
    audit["export"] = report
    audit["provenance"] = provenance
    write_audit(stem, audit)
    plt.close(fig)
    return report


def new_figure(height_px: int, **kwargs):
    return plt.figure(figsize=(WIDTH_IN, height_px / DPI), dpi=DPI, **kwargs)
