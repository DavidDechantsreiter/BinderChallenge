"""Render the Complexa pilot PDBs as annotated 3D views and animations."""

from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from biotite.structure import annotate_sse
from biotite.structure.io import pdbx
from biotite.structure.io.pdb import PDBFile
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter
from skimage.measure import marching_cubes


ROOT = Path(__file__).resolve().parent
HELIX = "#176b93"
STRAND = "#7451ab"
LOOP = "#89a7b9"
LIGAND = "#a67508"
CLASH = "#bf334d"
REFERENCE = (ROOT.parent.parent / "related_work/2026-microcystein/02_structures"
             / "microcystein/1LCM.cif")


@dataclass
class Atom:
    name: str
    residue: str
    chain: str
    element: str
    xyz: tuple[float, float, float]

    @property
    def label(self):
        return f"{self.chain}:{self.residue} {self.name}"


def reference_bonds():
    reference = pdbx.get_structure(pdbx.CIFFile.read(REFERENCE), model=1, include_bonds=True)
    heavy = reference[~np.isin(reference.element, ["H", "D"])]
    labels = [(f"{name}{res_id}", atom) for name, res_id, atom in
              zip(heavy.res_name, heavy.res_id, heavy.atom_name)]
    bonds = [(labels[i], labels[j]) for i, j, _ in heavy.bonds.as_array()]
    return [*bonds, (("FGA6", "C"), ("FGA6", "OXT"))]


def load_complex(path, reference_links):
    atoms = []
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and line[76:78].strip() not in {"H", "D"}:
            atoms.append(Atom(
                line[12:16].strip(), line[17:20].strip() + line[22:26].strip(),
                line[21], line[76:78].strip(), tuple(float(line[i:i + 8]) for i in (30, 38, 46)),
            ))
    ligand = [atom for atom in atoms if atom.chain == "A"]
    protein = [atom for atom in atoms if atom.chain == "B"]
    coords = np.array([atom.xyz for atom in ligand])
    backbone = np.array([atom.xyz for atom in protein if atom.name == "CA"])
    binder = PDBFile.read(path.with_name(f"{path.stem}_binder.pdb")).get_structure(model=1)
    if not np.allclose(backbone, binder.coord[binder.atom_name == "CA"]):
        raise ValueError(f"Binder coordinates differ from the complex in {path}")
    ribbon = protein_ribbon(backbone, annotate_sse(binder))
    distances = np.linalg.norm(coords[:, None] - np.array([atom.xyz for atom in protein]), axis=2)
    pairs = sorted(((distances[i, j], ligand[i], protein[j]) for i, j in zip(*np.where(distances < 2))),
                   key=lambda pair: pair[0])
    indices = {(atom.residue, atom.name): i for i, atom in enumerate(ligand)}
    ligand_bonds = [(indices[a], indices[b]) for a, b in reference_links]
    return ligand, coords, backbone, ligand_bonds, ribbon, molecular_volume(protein), pairs


def molecular_volume(atoms, spacing=0.55):
    """Sample the union of heavy-atom van der Waals spheres on a grid."""
    radii = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}
    coords = np.array([atom.xyz for atom in atoms])
    origin = coords.min(axis=0) - 3
    shape = np.ceil((coords.max(axis=0) + 3 - origin) / spacing).astype(int) + 1
    occupied = np.zeros(shape, dtype=bool)
    for atom, center in zip(atoms, (coords - origin) / spacing):
        radius = radii[atom.element] / spacing
        low = np.floor(center - radius).astype(int)
        high = np.ceil(center + radius).astype(int) + 1
        indices = np.ogrid[low[0]:high[0], low[1]:high[1], low[2]:high[2]]
        inside = sum((index - value) ** 2 for index, value in zip(indices, center)) <= radius**2
        occupied[tuple(slice(a, b) for a, b in zip(low, high))] |= inside
    vertices, faces, _, _ = marching_cubes(
        gaussian_filter(occupied.astype(np.float32), 0.9),
        level=0.45, spacing=(spacing,) * 3, step_size=2,
    )
    return float(occupied.sum() * spacing**3), vertices + origin, faces


def protein_ribbon(backbone, secondary):
    """Sweep an elliptical ribbon along a spline through the Cα positions."""
    if len(backbone) != len(secondary):
        raise ValueError("Secondary structure must have one assignment per Cα.")
    sample = np.linspace(0, len(backbone) - 1, 5 * (len(backbone) - 1) + 1)
    spline = CubicSpline(np.arange(len(backbone)), backbone, axis=0)
    centers = spline(sample)
    tangents = spline(sample, 1)
    tangents /= np.linalg.norm(tangents, axis=1, keepdims=True)
    residue = np.rint(sample).astype(int)
    colors = {"a": HELIX, "b": STRAND, "c": LOOP, "": LOOP}
    widths = {"a": 1.05, "b": 1.35, "c": 0.27, "": 0.27}
    width = np.interp(sample, np.arange(len(backbone)), [widths[code] for code in secondary])
    sides = 8
    angles = np.arange(sides) * 2 * np.pi / sides
    normal = np.eye(3)[np.argmin(np.abs(tangents[0]))]
    rings = []
    for center, tangent, radius in zip(centers, tangents, width):
        normal -= tangent * np.dot(normal, tangent)
        normal /= np.linalg.norm(normal)
        binormal = np.cross(tangent, normal)
        rings.append(center + radius * np.cos(angles)[:, None] * normal
                     + 0.23 * np.sin(angles)[:, None] * binormal)

    faces, face_colors = [], []
    for step, code in enumerate(secondary[residue[:-1]]):
        for j in range(sides):
            a = step * sides + j
            b = step * sides + (j + 1) % sides
            faces.extend([(a, b, a + sides), (b, b + sides, a + sides)])
            face_colors.extend([colors[code], colors[code]])
    return np.concatenate(rings), np.array(faces), face_colors


def draw_ribbon(ax, ribbon, opacity=1):
    vertices, faces, colors = ribbon
    ax.add_collection3d(Poly3DCollection(vertices[faces], facecolors=colors,
                                         linewidth=0, shade=True, alpha=opacity))


def approach_offset(coords, backbone):
    origin = coords.mean(axis=0)
    displacement = origin - backbone.mean(axis=0)
    return displacement * (17 / np.linalg.norm(displacement))


def pocket_view(coords, backbone, pairs):
    """Choose a view that exposes the ligand and the closest contact markers."""
    markers = np.array([np.mean([a.xyz, b.xyz], axis=0) for _, a, b in pairs])
    targets = np.vstack([coords, markers])
    offsets = backbone[None, :, :] - targets[:, None, :]
    lengths = np.sum(offsets**2, axis=2)
    outward = approach_offset(coords, backbone) / 17
    best = (-np.inf, 0, 0, outward)
    for elevation in range(-45, 46, 15):
        pitch = np.deg2rad(elevation)
        for azimuth in range(-180, 180, 10):
            yaw = np.deg2rad(azimuth)
            eye = np.array([np.cos(pitch) * np.cos(yaw), np.cos(pitch) * np.sin(yaw), np.sin(pitch)])
            depth = offsets @ eye
            hidden = np.any((depth > 0.5) & (lengths - depth**2 < 2.1**2), axis=1)
            projected = (coords - coords.mean(axis=0)) - np.outer((coords - coords.mean(axis=0)) @ eye, eye)
            spread = np.linalg.svd(projected, compute_uv=False)[1]
            score = (1 - hidden[:len(coords)].mean()) + 1.1 * (1 - hidden[len(coords):].mean())
            score += 0.015 * spread + 0.16 * (eye @ outward)
            if score > best[0]:
                best = (score, elevation, azimuth, eye)
    return best[1:]


def render_summary(models, closeup=False):
    figure = plt.figure(figsize=(10, 4.75), facecolor="#f5f8fa")
    for number, (_, coords, backbone, bonds, ribbon, _, pairs) in enumerate(models, 1):
        ax = figure.add_subplot(1, len(models), number, projection="3d", facecolor="#f5f8fa")
        if closeup:
            vertices, faces, colors = ribbon
            centers = vertices[faces].mean(axis=1)
            nearby = np.min(np.linalg.norm(centers[:, None] - coords, axis=-1), axis=1) < 5.5
            draw_ribbon(ax, (vertices, faces[nearby], np.array(colors)[nearby]), opacity=0.30)
            extent = np.vstack([coords, [atom.xyz for _, _, atom in pairs]])
        else:
            draw_ribbon(ax, ribbon)
            extent = np.vstack([ribbon[0], coords])
        for i, j in bonds:
            ax.plot(*coords[[i, j]].T, color=LIGAND, lw=2.2 if closeup else 1.6)
        ax.scatter(*coords.T, s=14 if closeup else 8, c=LIGAND, depthshade=False)
        for _, a, b in pairs:
            ax.plot(*np.array([a.xyz, b.xyz]).T, color=CLASH, lw=2.0,
                    marker="o" if closeup else None, markersize=3)
        distance, a, b = pairs[0]
        ax.scatter(*np.array([a.xyz, b.xyz]).T, s=30, c=CLASH, depthshade=False)
        ax.set_title(f"Design {number}  ·  {b.residue} {b.name} ↔ {a.residue} {a.name}  ·  {distance:.2f} Å",
                     fontsize=10, pad=-5)
        padding = 2.5 if closeup else 0.7
        for setter, values in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), extent.T):
            setter(values.min() - padding, values.max() + padding)
        ax.set_box_aspect(np.ptp(extent, axis=0) + 2 * padding, zoom=1.42)
        elevation, azimuth, _ = pocket_view(coords, backbone, pairs)
        ax.set_proj_type("ortho")
        ax.view_init(elev=elevation, azim=azimuth)
        ax.set_axis_off()
    title = "Binding-site close-ups" if closeup else "Pilot complexes"
    figure.text(0.5, 0.93, title, ha="center", fontsize=20, weight="bold", color="#203746")
    legend = [(0.15, "α helix", HELIX), (0.30, "β strand", STRAND), (0.43, "Loop", LOOP),
              (0.60, "Microcystin-LR", LIGAND), (0.81, "Pairs under 2 Å", CLASH)]
    for x, label, color in legend:
        figure.text(x, 0.065, label, ha="center", fontsize=10, weight="bold", color=color)
    figure.subplots_adjust(top=0.82, bottom=0.13, left=0.025, right=0.975, wspace=0.06)
    figure.savefig(ROOT / ("pocket-closeups.png" if closeup else "complexes.png"), dpi=155)
    plt.close(figure)


def render_animation(model):
    _, coords, backbone, bonds, ribbon, _, pairs = model
    displacement = approach_offset(coords, backbone)
    elevation, pocket_azimuth, _ = pocket_view(coords, backbone, pairs)
    extent = np.vstack([backbone, coords, coords + displacement])
    figure = plt.figure(figsize=(6.8, 4.8), facecolor="#f5f8fa")
    ax = figure.add_subplot(projection="3d", facecolor="#f5f8fa")
    ax.set_axis_off()
    ax.set_proj_type("ortho")
    ax.set_box_aspect(np.ptp(extent, axis=0), zoom=1.30)
    for setter, values in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), extent.T):
        setter(min(values) - 2, max(values) + 2)
    draw_ribbon(ax, ribbon)
    links = [ax.plot([], [], [], color=LIGAND, lw=1.8)[0] for _ in bonds]
    clash_marks = [
        ax.plot(*np.array([a.xyz, b.xyz]).T, color=CLASH, lw=2, marker="o", markersize=3)[0]
        for _, a, b in pairs
    ]
    figure.text(0.5, 0.94, "Design 2 · Illustrative approach", ha="center",
                fontsize=16, weight="bold", color="#203746")
    caption = figure.text(0.5, 0.08, "", ha="center", fontsize=11, color="#203746")

    def frame(state):
        offset, azimuth = state
        ax.view_init(elev=elevation + 5 * np.sin(np.deg2rad(azimuth - pocket_azimuth)), azim=azimuth)
        moved = coords + offset * displacement
        for artist, (i, j) in zip(links, bonds):
            artist.set_data_3d(*moved[[i, j]].T)
        for mark in clash_marks:
            mark.set_visible(offset == 0)
        caption.set_text("Final pose · 5 pairs under 2 Å" if offset == 0 else
                         "Binder in blue · ligand in gold")
        return [*links, *clash_marks, caption]

    progress = np.linspace(0, 1, 60)
    approach = np.clip((progress - 0.08) / 0.57, 0, 1)
    offsets = 1 - approach**2 * (3 - 2 * approach)
    states = list(zip(offsets, pocket_azimuth - 42 + 65 * progress))
    animation = FuncAnimation(figure, frame, frames=states, interval=85, repeat=True)
    animation.save(ROOT / "binding-design-2.gif", writer=PillowWriter(fps=12), dpi=95)
    plt.close(figure)


def render_tertiary(models):
    figure = plt.figure(figsize=(9, 4.8), facecolor="#f5f8fa")
    figure.text(0.5, 0.94, "Tertiary folds and relative molecular volume", ha="center",
                fontsize=17, weight="bold", color="#203746")
    figure.text(0.5, 0.06, "Same Å scale  ·  ribbons show the fold  ·  translucent shells show occupied volume",
                ha="center", fontsize=10, color="#203746")
    radius = max(np.max(np.linalg.norm(model[-2][1] - model[2].mean(axis=0), axis=1))
                 for model in models) + 1
    baseline = models[0][-2][0]
    axes = []
    for number, (_, coords, backbone, _, ribbon, volume, pairs) in enumerate(models, 1):
        amount, surface, triangles = volume
        ax = figure.add_subplot(1, 2, number, projection="3d", facecolor="#f5f8fa")
        center = backbone.mean(axis=0)
        for setter, value in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), center):
            setter(value - radius, value + radius)
        ax.set_box_aspect((1, 1, 1), zoom=1.24)
        ax.set_axis_off()
        draw_ribbon(ax, ribbon)
        ax.add_collection3d(Poly3DCollection(surface[triangles], facecolors="#7199ac",
                                             alpha=0.14, linewidth=0))
        ax.set_title(f"Design {number}  ·  {amount:,.0f} Å³  ·  {amount / baseline:.2f}× design 1",
                     fontsize=11)
        elevation, azimuth, _ = pocket_view(coords, backbone, pairs)
        axes.append((ax, elevation, azimuth))
        print(f"Design {number} heavy-atom volume: {amount:.1f} Å³ ({amount / baseline:.3f}×)")
    figure.subplots_adjust(left=0.01, right=0.99, top=0.85, bottom=0.12, wspace=0.01)

    def frame(rotation):
        for ax, elevation, azimuth in axes:
            ax.view_init(elev=elevation, azim=azimuth + rotation)

    frame(0)
    figure.savefig(ROOT / "tertiary-volume.png", dpi=145)
    animation = FuncAnimation(figure, frame, frames=np.linspace(0, 360, 60, endpoint=False))
    animation.save(ROOT / "tertiary-volume.gif", writer=PillowWriter(fps=12), dpi=85)
    plt.close(figure)


def main():
    paths = [path for path in sorted((ROOT / "structures").glob("job_*/*.pdb"))
             if not path.stem.endswith("_binder")]
    if len(paths) != 2:
        raise SystemExit(f"Expected two complex PDBs, found {len(paths)}.")
    reference_links = reference_bonds()
    models = [load_complex(path, reference_links) for path in paths]
    for number, model in enumerate(models, 1):
        for distance, ligand, protein in model[-1]:
            print(f"{number}: {protein.label} ↔ {ligand.label} {distance:.3f} Å")
    render_summary(models)
    render_summary(models, closeup=True)
    render_animation(models[1])
    render_tertiary(models)


if __name__ == "__main__":
    main()
