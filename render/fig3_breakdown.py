"""Figure 3: share of comparisons within the repeat range by instrument,
recorded collision energy, retained fragment count and precursor m/z.

Output: figures/overlap-by-measurement.png and .svg (3600 x 2880 px).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import common as C

from matplotlib.lines import Line2D

HEIGHT_PX = 2880
PANELS = [
    ("instrument", "Instrument"),
    ("collision_energy", "Recorded collision energy"),
    ("fragment_peaks", "Retained fragment peaks"),
    ("precursor_mz", "Precursor m/z (Th)"),
]
XLABEL = "Comparisons within repeat range (%)"
DOT_PT = 8.0


def independent_counts(iso: pd.DataFrame, threshold: float) -> dict:
    """Recount each category from the unit table to cross-check breakdowns.csv."""
    within = iso["score"] >= threshold
    e = iso["collision_energy"]
    k = iso["fragment_peak_count"]
    m = iso["precursor_mz"]
    masks = {
        ("instrument", "Orbitrap"): iso["instrument_type"] == "Orbitrap",
        ("instrument", "QTOF"): iso["instrument_type"] == "QTOF",
        ("collision_energy", "0"): e == 0,
        ("collision_energy", "(0,20]"): (e > 0) & (e <= 20),
        ("collision_energy", "(20,30]"): (e > 20) & (e <= 30),
        ("collision_energy", "(30,40]"): (e > 30) & (e <= 40),
        ("collision_energy", "(40,60]"): (e > 40) & (e <= 60),
        ("collision_energy", "over 60"): e > 60,
        ("fragment_peaks", "0 to 2"): k <= 2,
        ("fragment_peaks", "3 to 5"): (k >= 3) & (k <= 5),
        ("fragment_peaks", "6 to 10"): (k >= 6) & (k <= 10),
        ("fragment_peaks", "11 to 20"): (k >= 11) & (k <= 20),
        ("fragment_peaks", "21 to 50"): (k >= 21) & (k <= 50),
        ("fragment_peaks", "over 50"): k > 50,
        ("precursor_mz", "up to 200"): m <= 200,
        ("precursor_mz", "(200,300]"): (m > 200) & (m <= 300),
        ("precursor_mz", "(300,400]"): (m > 300) & (m <= 400),
        ("precursor_mz", "(400,500]"): (m > 400) & (m <= 500),
        ("precursor_mz", "over 500"): m > 500,
    }
    return {key: (int(mask.sum()), int((mask & within).sum())) for key, mask in masks.items()}


def main(formats: list[str]) -> None:
    br = C.read("breakdowns.csv")
    thr = C.read("threshold.csv")
    iso = C.read("isomer_scores.csv")
    threshold = float(thr["threshold"].iloc[0])
    overall = float(thr["within_percent"].iloc[0])
    C.check_scores(iso["score"], "isomer scores")

    br = br.sort_values(["dimension", "order"], kind="stable")
    assert len(br) == 19
    recomputed = 100.0 * br["within"] / br["n"]
    assert np.allclose(br["percent"], recomputed, rtol=0, atol=1e-9)
    assert br["percent"].max() < 50
    assert abs(overall - 100.0 * thr["within_units"].iloc[0] / thr["isomer_units"].iloc[0]) < 1e-9
    recount = independent_counts(iso, threshold)
    mismatches = []
    for row in br.itertuples():
        got = recount.get((row.dimension, row.label))
        if got != (row.n, row.within):
            mismatches.append([row.dimension, row.label, row.n, row.within, got])
    if mismatches:
        raise ValueError(f"Breakdown counts disagree with the unit table: {mismatches}")
    for dim in br["dimension"].unique():
        assert br.loc[br["dimension"] == dim, "n"].sum() == len(iso)

    labels = {"Overall"} | {p[1] for p in PANELS} | {XLABEL} | set(br["label"])

    with C.plt.rc_context(C.RC):
        fig = C.new_figure(HEIGHT_PX, layout="constrained")
        fig.get_layout_engine().set(w_pad=0.06, h_pad=0.08, wspace=0.06, hspace=0.06)
        sizes = [int((br["dimension"] == d).sum()) for d, _ in PANELS]
        top = max(sizes[0], sizes[1])
        bottom = max(sizes[2], sizes[3])
        grid = fig.add_gridspec(2, 2, height_ratios=[top + 0.6, bottom + 0.6])
        axes = [fig.add_subplot(grid[i // 2, i % 2]) for i in range(4)]
        for ax in axes[1:]:
            ax.sharex(axes[0])

        drawn = []
        for ax, (dim, name), rows_in_slot in zip(axes, PANELS, [top, top, bottom, bottom]):
            rows = br[br["dimension"] == dim].sort_values("order", kind="stable")
            ypos = np.arange(len(rows))
            ax.set_yticks(ypos, rows["label"].tolist())
            for yv in ypos:
                ax.axhline(yv, color=C.GRID, lw=0.6, zorder=0)
            ax.axvline(overall, color=C.GREY, lw=1.3, ls=(0, (5, 3)), zorder=1)
            ax.scatter(
                rows["percent"],
                ypos,
                s=DOT_PT**2,
                marker="o",
                facecolor=C.BLACK,
                edgecolor="white",
                linewidths=0.6,
                zorder=3,
                clip_on=False,
            )
            drawn.extend(zip([dim] * len(rows), rows["label"], rows["percent"]))
            # Same vertical spacing per category in every panel.
            pad = 0.6
            ax.set_ylim(len(rows) - 1 + pad + (rows_in_slot - len(rows)) / 2,
                        -pad - (rows_in_slot - len(rows)) / 2)
            ax.set_ylabel(name)
            ax.tick_params(axis="y", length=0, pad=6)
            ax.spines["left"].set_visible(False)
        axes[0].set_xlim(0, 50)
        axes[0].set_xticks(np.arange(0, 51, 10))
        for ax in axes[:2]:
            ax.tick_params(axis="x", labelbottom=True)
        fig.supxlabel(XLABEL, fontsize=C.RC["axes.labelsize"])

        # The instrument panel's lower rows are empty: the key covers no data.
        axes[0].legend(
            [Line2D([], [], color=C.GREY, lw=1.3, ls=(0, (5, 3)))],
            ["Overall"],
            loc="lower right",
            handlelength=2.6,
            borderaxespad=0.4,
        )

        fig.canvas.draw()
        audit = {
            "checks": {
                "rows": len(br),
                "percent_equals_100_within_over_n": True,
                "overall_percent": overall,
                "recount_mismatches_against_isomer_scores": mismatches,
            },
            "drawn_points": [[d, lab, float(p)] for d, lab, p in drawn],
        }
        audit["text"] = C.audit_text(fig, labels)
        assert not audit["text"]["disallowed"], audit["text"]["disallowed"]
        C.export(
            fig,
            "overlap-by-measurement",
            formats,
            provenance={
                "raw_data": [
                    "figures/data/breakdowns.csv",
                    "figures/data/threshold.csv",
                ],
                "transformations": ["none: percent column plotted as supplied"],
                "uncertainty": "none drawn",
                "missing_data": "none present",
            },
            audit=audit,
        )


if __name__ == "__main__":
    main(C.arguments(__doc__, ("png", "svg")))
