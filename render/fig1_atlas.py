"""Figure 1: every isomer comparison by precursor m/z and similarity, with the
empirical repeat-score survival curves on a shared similarity axis.

Output: figures/isomer-atlas.png (3600 x 2000 px).
"""

from __future__ import annotations

import numpy as np

import common as C

from matplotlib.collections import PathCollection
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from matplotlib.markers import MarkerStyle
from matplotlib.transforms import IdentityTransform

HEIGHT_PX = 2000
POINT_PT = 2.3  # circle diameter in points; triangles have the same area
POINT_ALPHA = 0.34
DIAMOND_PT = 3.6
FOLDS = [
    ("train", "Training", "-"),
    ("val", "Validation", (0, (4.5, 2.2))),
    ("test", "Test", (0, (0.1, 2.0))),
]
ALLOWED = {
    "Precursor m/z (Th)",
    "Fragment similarity",
    "Repeats at or above (%)",
    "Orbitrap",
    "QTOF",
    "Equal lists (tolerance)",
    "Training",
    "Validation",
    "Test",
}


def marker_path(symbol: str, scale: float = 1.0):
    style = MarkerStyle(symbol)
    return style.get_path().transformed(style.get_transform().scale(scale))


def main(formats: list[str]) -> None:
    iso = C.read("isomer_scores.csv")
    rep = C.read("repeat_scores.csv")
    thr = C.read("threshold.csv")
    threshold = float(thr["threshold"].iloc[0])

    C.check_scores(iso["score"], "isomer scores")
    C.check_scores(rep["score"], "repeat scores")
    n_units = len(iso)
    checks = {
        "units": n_units,
        "units_expected": int(thr["isomer_units"].iloc[0]),
        "within": int((iso["score"] >= threshold).sum()),
        "within_expected": int(thr["within_units"].iloc[0]),
        "tolerance_flagged": int(iso["identical_peak_list"].sum()),
        "unit_key_unique": bool(not iso.duplicated(["query_id", "reference_molecule"]).any()),
        "training_repeats": int((rep["fold"] == "train").sum()),
        "training_repeats_expected": int(thr["training_repeats"].iloc[0]),
        "precursor_mz_range": [float(iso["precursor_mz"].min()), float(iso["precursor_mz"].max())],
        "score_range": [float(iso["score"].min()), float(iso["score"].max())],
        "instruments": iso["instrument_type"].value_counts().to_dict(),
        "fold_sizes": rep["fold"].value_counts().to_dict(),
    }
    assert checks["units"] == checks["units_expected"] == 62535
    assert checks["within"] == checks["within_expected"]
    assert checks["tolerance_flagged"] == 897
    assert checks["unit_key_unique"]
    assert checks["training_repeats"] == checks["training_repeats_expected"]
    assert 0 <= checks["precursor_mz_range"][0] and checks["precursor_mz_range"][1] <= 1000
    assert set(checks["instruments"]) == {"Orbitrap", "QTOF"}

    units = iso.sort_values(["query_id", "reference_molecule"], kind="stable")
    x = units["precursor_mz"].to_numpy(dtype=float)
    y = np.clip(units["score"].to_numpy(dtype=float), 0.0, 1.0)  # drawing only
    orbitrap = (units["instrument_type"] == "Orbitrap").to_numpy()
    flagged = units["identical_peak_list"].to_numpy(dtype=bool)

    with C.plt.rc_context(C.RC):
        fig = C.new_figure(HEIGHT_PX, layout="constrained")
        fig.get_layout_engine().set(w_pad=0.06, h_pad=0.06, wspace=0.02)
        grid = fig.add_gridspec(1, 2, width_ratios=[3.05, 1.0])
        ax = fig.add_subplot(grid[0, 0])
        axr = fig.add_subplot(grid[0, 1], sharey=ax)

        for a in (ax, axr):
            a.axhspan(threshold, 1.0, color=C.GREEN, alpha=0.07, lw=0, zorder=0)
            a.grid(True, axis="y", zorder=0)
        ax.grid(True, axis="x", zorder=0)

        # One collection in unit-key order so neither instrument is painted last.
        circle = marker_path("o")
        triangle = marker_path("^", np.sqrt((np.pi / 4) / 0.5))  # equal area
        paths = [circle if o else triangle for o in orbitrap]
        colours = np.where(
            orbitrap[:, None],
            np.array(to_rgba(C.BLUE, POINT_ALPHA))[None, :],
            np.array(to_rgba(C.ORANGE, POINT_ALPHA))[None, :],
        )
        points = PathCollection(
            paths,
            sizes=[POINT_PT**2],
            offsets=np.column_stack([x, y]),
            offset_transform=ax.transData,
            facecolors=colours,
            edgecolors="none",
            linewidths=0.0,
            zorder=2,
        )
        points.set_transform(IdentityTransform())  # marker sizes in points, as scatter does
        ax.add_collection(points, autolim=False)
        ax.scatter(
            x[flagged],
            y[flagged],
            marker="D",
            s=DIAMOND_PT**2,
            facecolors="none",
            edgecolors=C.CHARCOAL,
            linewidths=0.5,
            zorder=3,
        )
        for a in (ax, axr):
            a.axhline(threshold, color=C.BLACK, lw=1.0, ls=(0, (5, 3)), zorder=4)

        curve_checks = {}
        for fold, label, style in FOLDS:
            scores = rep.loc[rep["fold"] == fold, "score"].to_numpy(dtype=float)
            t, pct = C.survival_vertices(scores, scores.size, 0.0, 1.0)
            axr.plot(
                pct,
                np.clip(t, 0.0, 1.0),
                color=C.GREEN,
                lw=1.5,
                ls=style,
                dash_capstyle="round",
                solid_joinstyle="miter",
                clip_on=False,
                zorder=3,
                label=label,
            )
            at_cut = C.evaluate_survival(t, pct, threshold)
            exact = 100.0 * np.count_nonzero(scores >= threshold) / scores.size
            assert abs(at_cut - exact) < 1e-12
            curve_checks[fold] = {
                "n": int(scores.size),
                "vertices": int(t.size),
                "percent_at_or_above_cutoff": exact,
            }

        ax.set_xlim(0, 1000)
        ax.set_ylim(-0.025, 1.025)
        ax.set_xticks(np.arange(0, 1001, 200))
        ax.set_yticks(np.linspace(0, 1, 6))
        ax.set_xlabel("Precursor m/z (Th)")
        ax.set_ylabel("Fragment similarity")
        axr.set_xlim(0, 100)
        axr.set_xticks([0, 50, 100])
        axr.set_xlabel("Repeats at or\nabove (%)")
        axr.tick_params(axis="y", labelleft=False)

        key = [
            Line2D([], [], ls="none", marker="o", ms=7, mfc=C.BLUE, mec="none"),
            Line2D([], [], ls="none", marker="^", ms=7.5, mfc=C.ORANGE, mec="none"),
            Line2D([], [], ls="none", marker="D", ms=6.5, mfc="none", mec=C.CHARCOAL, mew=1.1),
        ]
        ax.legend(
            key,
            ["Orbitrap", "QTOF", "Equal lists (tolerance)"],
            loc="lower left",
            bbox_to_anchor=(-0.01, 1.0),
            ncol=3,
            handletextpad=0.3,
            columnspacing=1.5,
            borderaxespad=0.15,
        )
        axr.legend(
            loc="lower left",
            bbox_to_anchor=(0.0, 0.02),
            handlelength=2.6,
            handletextpad=0.5,
            borderaxespad=0.2,
        )

        fig.canvas.draw()
        audit = {"checks": checks, "repeat_curves": curve_checks, "threshold": threshold}
        audit["text"] = C.audit_text(fig, ALLOWED)
        assert not audit["text"]["disallowed"], audit["text"]["disallowed"]
        C.export(
            fig,
            "isomer-atlas",
            formats,
            provenance={
                "raw_data": [
                    "figures/data/isomer_scores.csv",
                    "figures/data/repeat_scores.csv",
                    "figures/data/threshold.csv",
                ],
                "transformations": [
                    "all 62,535 units drawn, ordered by (query_id, reference_molecule)",
                    "drawing y clipped to [0, 1]; stored scores unchanged",
                    "repeat curves: exact tie-aware empirical survival per fold",
                ],
                "uncertainty": "none drawn",
                "missing_data": "none present",
            },
            audit=audit,
        )


if __name__ == "__main__":
    main(C.arguments(__doc__, ("png",)))
