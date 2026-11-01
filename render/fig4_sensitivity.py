"""Figure 4: exact empirical share of comparisons at or above each similarity
cutoff, for all comparisons and with tolerance-equal peak lists excluded
(original denominator kept), with the four checked repeat-quantile cutoffs.

Output: figures/threshold-sensitivity.png (3600 x 2100 px).
"""

from __future__ import annotations

import numpy as np

import common as C

HEIGHT_PX = 2100
ALLOWED = {
    "Similarity cutoff",
    "Comparisons at or above cutoff (%)",
    "Lower training-repeat quantile (%)",
    "All comparisons",
    "Equal lists excluded",
}


def main(formats: list[str]) -> None:
    iso = C.read("isomer_scores.csv")
    sens = C.read("sensitivity.csv").sort_values("quantile", kind="stable")
    thr = C.read("threshold.csv")
    main_cut = float(thr["threshold"].iloc[0])
    C.check_scores(iso["score"], "isomer scores")

    scores = iso["score"].to_numpy(dtype=float)
    equal = iso["identical_peak_list"].to_numpy(dtype=bool)
    n = scores.size
    t_all, p_all = C.survival_vertices(scores, n, 0.0, 1.0)
    t_kept, p_kept = C.survival_vertices(scores[~equal], n, 0.0, 1.0)

    point_checks = []
    for row in sens.itertuples():
        on_curve = C.evaluate_survival(t_all, p_all, row.threshold)
        count = int(np.count_nonzero(scores >= row.threshold))
        assert row.n == n and count == row.within
        assert abs(on_curve - row.percent) < 1e-9, (on_curve, row.percent)
        point_checks.append(
            {"quantile": row.quantile, "threshold": row.threshold, "curve": on_curve, "csv": row.percent}
        )
    assert main_cut in set(sens["threshold"])
    kept_at_main = C.evaluate_survival(t_kept, p_kept, main_cut)
    kept_count = int(np.count_nonzero((scores >= main_cut) & ~equal))
    assert kept_count == 10657 and f"{kept_at_main:.2f}" == "17.04", kept_at_main
    assert equal.sum() == 897

    tick_labels = [f"{100 * q:g}" for q in sens["quantile"]]
    assert tick_labels == ["0.5", "1", "2", "5"]
    allowed = ALLOWED | set(tick_labels)

    with C.plt.rc_context(C.RC):
        fig = C.new_figure(HEIGHT_PX, layout="constrained")
        fig.get_layout_engine().set(w_pad=0.06, h_pad=0.06)
        ax = fig.add_subplot()
        ax.grid(True, axis="y", zorder=0)
        for row in sens.itertuples():
            main_rule = row.threshold == main_cut
            ax.axvline(
                row.threshold,
                color=C.GREY if main_rule else C.LIGHT_GREY,
                lw=1.2 if main_rule else 0.8,
                zorder=1,
            )
        ax.plot(
            np.clip(t_all, 0.0, 1.0), p_all,
            color=C.BLUE, lw=1.8, solid_joinstyle="miter", zorder=3,
            clip_on=False, label="All comparisons",
        )
        ax.plot(
            np.clip(t_kept, 0.0, 1.0), p_kept,
            color=C.CHARCOAL, lw=1.8, ls=(0, (4, 2.2)), zorder=4,
            clip_on=False, label="Equal lists excluded",
        )
        ax.set_xlim(-0.01, 1.01)
        ax.set_ylim(-1.5, 101.5)
        ax.set_xticks(np.linspace(0, 1, 6))
        ax.set_yticks(np.arange(0, 101, 20))
        ax.set_xlabel("Similarity cutoff")
        ax.set_ylabel("Comparisons at or above cutoff (%)")
        top = ax.secondary_xaxis("top")
        top.set_xticks(sens["threshold"].to_numpy(), tick_labels)
        top.set_xlabel("Lower training-repeat quantile (%)", loc="right")
        top.tick_params(axis="x", pad=3)
        ax.spines["top"].set_visible(False)
        ax.legend(
            loc="upper left",
            bbox_to_anchor=(0.36, 1.0),
            handlelength=2.8,
            borderaxespad=0.3,
        )

        fig.canvas.draw()
        audit = {
            "checks": {
                "units": n,
                "tolerance_flagged": int(equal.sum()),
                "sensitivity_points_on_first_curve": point_checks,
                "equal_lists_excluded_at_main_cutoff": {"count": kept_count, "percent": kept_at_main},
                "vertices": {"all": int(t_all.size), "kept": int(t_kept.size)},
            }
        }
        audit["text"] = C.audit_text(fig, allowed)
        assert not audit["text"]["disallowed"], audit["text"]["disallowed"]
        C.export(
            fig,
            "threshold-sensitivity",
            formats,
            provenance={
                "raw_data": [
                    "figures/data/isomer_scores.csv",
                    "figures/data/sensitivity.csv",
                    "figures/data/threshold.csv",
                ],
                "transformations": [
                    "exact tie-aware empirical survival over all 62,535 units",
                    "second curve drops identical_peak_list units, denominator stays 62,535",
                    "drawing x clipped to [0, 1]; stored scores unchanged",
                ],
                "uncertainty": "none drawn",
                "missing_data": "none present",
            },
            audit=audit,
        )


if __name__ == "__main__":
    # PNG preserves every step; SVG paths of this size can render incorrectly.
    main(C.arguments(__doc__, ("png",)))
