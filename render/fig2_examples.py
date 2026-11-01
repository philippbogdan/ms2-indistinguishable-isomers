"""Figure 2: the three supplied example pairs, each as two RDKit structures
(query above reference) beside a mirror spectrum of retained peaks.

Output: figures/isomer-examples.png and .svg.

Structures come from the supplied SMILES through RDKit only. The reference
layout is generated to match the query layout on their maximum common
substructure (coordinates only; connectivity is checked unchanged). All six
molecules use one fixed bond length. The PNG composites RDKit's Cairo
renders pixel for pixel; the SVG nests RDKit's own SVG output.
"""

from __future__ import annotations

import io
import math
import re

import common as C

import numpy as np
from matplotlib.ticker import FuncFormatter, MultipleLocator
from PIL import Image
from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor, rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D

PAIR_ORDER = ["triazines", "sulfones", "flavonoids"]
BOND_PX = 120  # native px per bond, identical for all six molecules
ATOM_FONT_PX = 76  # native px; 16.9 px em, about 12.3 px cap height at 800 px
BOND_LINE_PX = 6
NAME_GAP_PX = 18
NAME_BAND_PX = 110
ROW_GAP_PX = 260
TOP_PX = 30
BOTTOM_PX = 250
LEFT_PX = 30
COLUMN_PX = 1320
AXES_LEFT_PX = 1700
RIGHT_PX = 90
STICK_LW = 1.5
YLIM = 106.0


def round_up(value: float, step: float = 25.0) -> float:
    return math.ceil(value / step) * step


def depict_pair(query_smiles: str, reference_smiles: str):
    query = Chem.MolFromSmiles(query_smiles)
    reference_input = Chem.MolFromSmiles(reference_smiles)
    rdDepictor.SetPreferCoordGen(False)  # RDKit's depictor: honours canonOrient and exact templates
    rdDepictor.Compute2DCoords(query, canonOrient=True)
    mcs = rdFMCS.FindMCS(
        [query, reference_input],
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrder,
        timeout=5,
    )
    if mcs.canceled:
        raise ValueError("Shared-substructure search timed out")
    core = Chem.MolFromSmarts(mcs.smartsString)
    q_matches = query.GetSubstructMatches(core, uniquify=False, maxMatches=200)
    r_matches = reference_input.GetSubstructMatches(core, uniquify=False, maxMatches=200)

    def matched_offset(mol, q_match, r_match) -> float:
        q_conf, r_conf = query.GetConformer(), mol.GetConformer()
        return max((q_conf.GetAtomPosition(a) - r_conf.GetAtomPosition(b)).Length() for a, b in zip(q_match, r_match))

    # Symmetric cores have several valid atom correspondences; keep the one
    # whose overlay is most compact (exact template placement first).
    best = None
    candidates = 0
    for q_match in q_matches:
        for r_match in r_matches:
            trial = Chem.Mol(reference_input)
            # (mol, reference, atomMap of (reference idx, mol idx), confId, forceRDKit)
            rdDepictor.GenerateDepictionMatching2DStructure(
                trial, query, list(zip(q_match, r_match)), -1, True
            )
            union = np.vstack([coords(query), coords(trial)])
            area = float(np.ptp(union[:, 0]) * np.ptp(union[:, 1]))
            offset = matched_offset(trial, q_match, r_match)
            key = (offset > 1e-3, area)
            candidates += 1
            if best is None or key < best[0]:
                best = (key, trial, q_match, r_match, offset)
    aligned = best is not None
    if aligned:
        _, reference, q_match, r_match, deviation = best
    else:
        reference = Chem.Mol(reference_input)
        rdDepictor.Compute2DCoords(reference, canonOrient=True)
        q_match, r_match, deviation = (), (), float("nan")
    # Identical shared rings should show the same Kekule form. Kekulise both,
    # copy the query's bond orders onto the reference's matched bonds, and keep
    # that only if RDKit sanitises it to the very same molecule (canonical
    # SMILES); otherwise keep RDKit's own form. Drawing then skips
    # re-kekulisation so the chosen form is what appears.
    canonical = {
        "query": Chem.MolToSmiles(Chem.MolFromSmiles(query_smiles)),
        "reference": Chem.MolToSmiles(Chem.MolFromSmiles(reference_smiles)),
    }
    kq, kr = Chem.Mol(query), Chem.Mol(reference)
    Chem.Kekulize(kq, clearAromaticFlags=True)
    Chem.Kekulize(kr, clearAromaticFlags=True)
    pairs_map = dict(zip(q_match, r_match))
    harmonised = Chem.RWMol(kr)
    changed = 0
    for bond in kq.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in pairs_map and b in pairs_map and query.GetBondBetweenAtoms(a, b).GetIsAromatic():
            target = harmonised.GetBondBetweenAtoms(pairs_map[a], pairs_map[b])
            if target is not None and target.GetBondType() != bond.GetBondType():
                target.SetBondType(bond.GetBondType())
                changed += 1
    harmonised = harmonised.GetMol()
    check = Chem.Mol(harmonised)
    try:
        Chem.SanitizeMol(check)
        same = Chem.MolToSmiles(check) == canonical["reference"]
    except Exception:  # noqa: BLE001 - any sanitisation failure means keep RDKit's form
        same = False
    if same:
        kr = harmonised
    kekule = {"reference_bonds_reassigned": changed if same else 0, "harmonised": bool(same)}
    query = rdMolDraw2D.PrepareMolForDrawing(kq, kekulize=False)
    reference = rdMolDraw2D.PrepareMolForDrawing(kr, kekulize=False)
    for mol, role in ((query, "query"), (reference, "reference")):
        assert mol.GetNumAtoms() == Chem.MolFromSmiles(query_smiles if role == "query" else reference_smiles).GetNumAtoms()
        copy = Chem.Mol(mol)
        Chem.SanitizeMol(copy)
        assert Chem.MolToSmiles(copy) == canonical[role]
    info = {
        "mcs_smarts": mcs.smartsString,
        "mcs_atoms": mcs.numAtoms,
        "aligned_on_mcs": aligned,
        "correspondences_tried": candidates,
        "kekule": kekule,
        "max_matched_atom_offset": deviation,
        "query_match": list(q_match),
        "reference_match": list(r_match),
    }
    return query, reference, info


def coords(mol) -> np.ndarray:
    conf = mol.GetConformer()
    return np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y] for i in range(mol.GetNumAtoms())])


def mean_bond_length(mol) -> float:
    xy = coords(mol)
    return float(np.mean([np.linalg.norm(xy[b.GetBeginAtomIdx()] - xy[b.GetEndAtomIdx()]) for b in mol.GetBonds()]))


def configure(drawer) -> None:
    opts = drawer.drawOptions()
    opts.useBWAtomPalette()
    opts.fixedBondLength = BOND_PX / 1.5  # RDKit: px per coordinate unit; its bonds are 1.5 units
    opts.bondLineWidth = BOND_LINE_PX
    opts.scaleBondWidth = False
    opts.minFontSize = ATOM_FONT_PX
    opts.maxFontSize = ATOM_FONT_PX
    opts.padding = 0.0
    opts.clearBackground = False
    opts.addStereoAnnotation = False
    opts.addAtomIndices = False
    opts.prepareMolsBeforeDrawing = False  # molecules arrive kekulised and prepared


def render(mol, width: int, height: int, kind: str):
    drawer = (
        rdMolDraw2D.MolDraw2DCairo(width, height)
        if kind == "png"
        else rdMolDraw2D.MolDraw2DSVG(width, height)
    )
    configure(drawer)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    pixels = np.array([[drawer.GetDrawCoords(i).x, drawer.GetDrawCoords(i).y] for i in range(mol.GetNumAtoms())])
    return drawer.GetDrawingText(), pixels


def pixel_scale(mol, pixels: np.ndarray) -> float:
    xy = coords(mol)
    ratios = [
        np.linalg.norm(pixels[b.GetBeginAtomIdx()] - pixels[b.GetEndAtomIdx()])
        / np.linalg.norm(xy[b.GetBeginAtomIdx()] - xy[b.GetEndAtomIdx()])
        for b in mol.GetBonds()
    ]
    return float(np.median(ratios))


def main(formats: list[str]) -> None:
    pairs = C.read("example_pairs.csv").set_index("pair_id").loc[PAIR_ORDER].reset_index()
    spectra = C.read("example_spectra.csv")
    peaks = C.read("example_peaks.csv")
    by_id = spectra.set_index("spectrum_id")

    checks = []
    data = []
    for pair in pairs.itertuples():
        q, r = by_id.loc[pair.query_id], by_id.loc[pair.reference_id]
        assert q["pair_id"] == r["pair_id"] == pair.pair_id
        assert q["role"] == "query" and r["role"] == "reference"
        assert q["name"] == pair.query_name and r["name"] == pair.reference_name
        assert q["molecule_key"] != r["molecule_key"]
        assert not pair.identical_peak_list and not pair.exact_peak_list
        sides = {}
        for role, row, spectrum_id, expected in (
            ("query", q, pair.query_id, pair.query_fragment_peaks),
            ("reference", r, pair.reference_id, pair.reference_fragment_peaks),
        ):
            rows = peaks[peaks["spectrum_id"] == spectrum_id].sort_values("peak_index", kind="stable")
            kept = rows[rows["retained"]]
            assert len(kept) == expected and len(kept) >= 8
            height = 100.0 * kept["relative_intensity"].to_numpy() / kept["relative_intensity"].max()
            sides[role] = {
                "name": row["name"],
                "smiles": row["smiles"],
                "precursor_mz": float(row["precursor_mz"]),
                "mz": kept["mz"].to_numpy(dtype=float),
                "height": height,
                "all_peaks": len(rows),
            }
        xmax = round_up(max(sides["query"]["precursor_mz"], sides["reference"]["precursor_mz"]))
        for side in sides.values():
            assert side["mz"].min() > 0 and side["mz"].max() <= xmax
        data.append({"pair": pair, "sides": sides, "xmax": xmax})
        checks.append(
            {
                "pair_id": pair.pair_id,
                "query": [pair.query_id, sides["query"]["name"], int(len(sides["query"]["mz"])), sides["query"]["all_peaks"]],
                "reference": [pair.reference_id, sides["reference"]["name"], int(len(sides["reference"]["mz"])), sides["reference"]["all_peaks"]],
                "formula": pair.formula,
                "condition": pair.condition,
                "score": pair.score,
                "x_limit": xmax,
            }
        )

    # Cross-check each pair against the unit table: same key, same score.
    iso = C.read("isomer_scores.csv")
    for entry, check in zip(data, checks):
        pair = entry["pair"]
        ref_key = by_id.loc[pair.reference_id, "molecule_key"]
        unit = iso[(iso["query_id"] == pair.query_id) & (iso["reference_molecule"] == ref_key)]
        assert len(unit) == 1 and unit["score"].iloc[0] == pair.score
        assert bool(unit["identical_peak_list"].iloc[0]) is False
        check["unit_table_score_matches"] = True

    # Structures.
    mols = []
    for entry in data:
        query, reference, info = depict_pair(entry["sides"]["query"]["smiles"], entry["sides"]["reference"]["smiles"])
        entry["mols"] = (query, reference)
        entry["mcs"] = info
        mols.extend([query, reference])

    # Probe each drawing on a generous canvas to measure the ink extent in px.
    probe_w, probe_h = 2400, 1600
    extents = []
    scales = []
    for mol in mols:
        png, pixels = render(mol, probe_w, probe_h, "png")
        alpha = np.asarray(Image.open(io.BytesIO(png)).convert("RGBA"))[..., 3]
        ys, xs = np.nonzero(alpha)
        extents.append((int(xs.max() - xs.min() + 1), int(ys.max() - ys.min() + 1)))
        scales.append(pixel_scale(mol, pixels))
    # One drawing scale for all six molecules: comparable bond lengths.
    assert max(scales) / min(scales) - 1 < 0.005, scales
    ink_w = max(e[0] for e in extents)
    ink_h = max(e[1] for e in extents)
    assert ink_w <= COLUMN_PX, (ink_w, COLUMN_PX)
    struct_h = int(ink_h + 40)
    canvas_w, canvas_h = COLUMN_PX + 400, struct_h + 200  # transparent margin for offsets

    row_h = 2 * (NAME_BAND_PX + struct_h)
    height_px = TOP_PX + 3 * row_h + 2 * ROW_GAP_PX + BOTTOM_PX
    width_px = int(C.WIDTH_IN * C.DPI)

    with C.plt.rc_context(C.RC):
        fig = C.new_figure(height_px)
        placements = []
        axes = []
        for index, entry in enumerate(data):
            row_top = height_px - TOP_PX - index * (row_h + ROW_GAP_PX)  # from bottom, px
            row_bottom = row_top - row_h
            ax = fig.add_axes(
                [
                    AXES_LEFT_PX / width_px,
                    row_bottom / height_px,
                    (width_px - AXES_LEFT_PX - RIGHT_PX) / width_px,
                    row_h / height_px,
                ]
            )
            axes.append(ax)
            sides = entry["sides"]
            ax.vlines(sides["query"]["mz"], 0, sides["query"]["height"], colors=C.BLUE, lw=STICK_LW, capstyle="butt", zorder=3)
            ax.vlines(sides["reference"]["mz"], 0, -sides["reference"]["height"], colors=C.ORANGE, lw=STICK_LW, capstyle="butt", zorder=3)
            ax.axhline(0, color="#222222", lw=0.8, zorder=2)
            ax.set_xlim(0, entry["xmax"])
            ax.set_ylim(-YLIM, YLIM)
            ax.xaxis.set_major_locator(MultipleLocator(50))
            ax.set_yticks([-100, -50, 0, 50, 100])
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{abs(v):g}"))
            ax.grid(True, axis="y", zorder=0)
            if index == 1:
                ax.set_ylabel("Relative intensity (%)")
            if index == len(data) - 1:
                ax.set_xlabel("m/z (Th)")

            # Structure boxes: query in the upper half of the row, reference in the lower half.
            centre_x = LEFT_PX + COLUMN_PX / 2
            query, reference = entry["mols"]
            q_box_top = row_top - NAME_BAND_PX
            r_box_bottom = row_bottom + NAME_BAND_PX
            boxes = {
                "query": (query, q_box_top - struct_h / 2),
                "reference": (reference, r_box_bottom + struct_h / 2),
            }
            # Common horizontal frame: matched scaffold atoms share x across the pair.
            union = np.vstack([coords(query), coords(reference)])
            frame_centre = 0.5 * (union[:, 0].min() + union[:, 0].max())
            for role, (mol, centre_y) in boxes.items():
                placements.append(
                    {"pair": entry["pair"].pair_id, "role": role, "mol": mol, "centre_x": centre_x, "centre_y": centre_y, "frame_centre": frame_centre}
                )
            fig.text(
                centre_x / width_px,
                (q_box_top + NAME_GAP_PX) / height_px,
                sides["query"]["name"],
                ha="center",
                va="bottom",
                fontsize=C.RC["font.size"],
            )
            fig.text(
                centre_x / width_px,
                (r_box_bottom - NAME_GAP_PX) / height_px,
                sides["reference"]["name"],
                ha="center",
                va="top",
                fontsize=C.RC["font.size"],
            )

        # Place each structure so that frame_centre maps to the column centre
        # and the structure's own vertical extent is centred in its box.
        images = []
        for item in placements:
            png, pixels = render(item["mol"], canvas_w, canvas_h, "png")
            rgba = np.asarray(Image.open(io.BytesIO(png)).convert("RGBA"))
            scale = pixel_scale(item["mol"], pixels)
            xy = coords(item["mol"])
            # canvas px of the frame centre (x) and of the ink centre (y)
            frame_px = pixels[0, 0] + (item["frame_centre"] - xy[0, 0]) * scale
            ys, xs = np.nonzero(rgba[..., 3])
            ink_mid_y = 0.5 * (ys.min() + ys.max())
            x0 = int(round(item["centre_x"] - frame_px))
            y_top = int(round(item["centre_y"] + ink_mid_y))  # figure px from bottom of canvas top edge
            y0 = y_top - canvas_h
            ink_left, ink_right = x0 + xs.min(), x0 + xs.max()
            assert ink_left >= 0 and ink_right <= AXES_LEFT_PX - 260, (item["pair"], item["role"], ink_left, ink_right)
            item.update({"x0": x0, "y0": y0, "scale_px_per_unit": scale, "ink_x": [int(ink_left), int(ink_right)]})
            images.append(fig.figimage(rgba, xo=x0, yo=y0, origin="upper", zorder=5))

        fig.canvas.draw()
        allowed = {"Relative intensity (%)", "m/z (Th)"} | set(spectra["name"])
        audit = {
            "rdkit": rdBase.rdkitVersion,
            "pairs": checks,
            "mcs": {entry["pair"].pair_id: entry["mcs"] for entry in data},
            "structures": [
                {k: v for k, v in item.items() if k != "mol"} | {"smiles_drawn": Chem.MolToSmiles(item["mol"])}
                for item in placements
            ],
            "bond_px": BOND_PX,
            "px_per_coordinate_unit": scales,
            "mean_bond_length_px": [s_ * mean_bond_length(m) for s_, m in zip(scales, mols)],
            "atom_font_px": ATOM_FONT_PX,
            "size_px": [width_px, height_px],
        }
        audit["text"] = C.audit_text(fig, allowed)
        assert not audit["text"]["disallowed"], audit["text"]["disallowed"]
        provenance = {
            "raw_data": [
                "figures/data/example_pairs.csv",
                "figures/data/example_spectra.csv",
                "figures/data/example_peaks.csv",
            ],
            "transformations": [
                "retained peaks only; each spectrum divided by its retained maximum, times 100",
                "reference spectrum drawn downwards; no peak merged, dropped or shifted",
                "structures: RDKit depiction of supplied SMILES, reference matched to query on MCS",
            ],
            "uncertainty": "none drawn",
            "missing_data": "none present",
        }
        if not audit["text"]["all_inside_canvas"]:
            raise ValueError("Molecular example labels extend outside the canvas")
        audit["provenance"] = provenance
        if "png" in formats:
            audit["export_png"] = C.save(fig, "isomer-examples", ["png"])

        if "svg" in formats:
            for image in images:
                image.remove()
            report_svg = C.save(fig, "isomer-examples", ["svg"])
            svg_path = C.FIGURES / "isomer-examples.svg"
            text = svg_path.read_text()
            px_to_pt = 72.0 / C.DPI
            nested = []
            for item in placements:
                svg, _ = render(item["mol"], canvas_w, canvas_h, "svg")
                opening = re.search(r"<svg\b[^>]*>", svg)
                inner = svg[opening.end() : svg.rfind("</svg>")]
                x_pt = item["x0"] * px_to_pt
                y_pt = (height_px - (item["y0"] + canvas_h)) * px_to_pt
                nested.append(
                    f'<svg x="{x_pt:.3f}" y="{y_pt:.3f}" width="{canvas_w * px_to_pt:.3f}" '
                    f'height="{canvas_h * px_to_pt:.3f}" viewBox="0 0 {canvas_w} {canvas_h}" '
                    f'overflow="visible" xmlns:rdkit="http://www.rdkit.org/xml">{inner}</svg>'
                )
            closing = text.rstrip().rfind("</svg>")
            text = text[:closing] + "\n".join(nested) + "\n</svg>\n"
            svg_path.write_text(text)
            report_svg["outputs"][0]["bytes"] = svg_path.stat().st_size
            audit["export_svg"] = report_svg

        C.write_audit("isomer-examples", audit)
        C.plt.close(fig)


if __name__ == "__main__":
    main(C.arguments(__doc__, ("png", "svg")))
