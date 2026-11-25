"""Summarize the interface geometry of a Proteina-Complexa pilot."""

import argparse
import csv
from pathlib import Path
import sys

import numpy as np


def read_heavy_atoms(path):
    atoms = {"A": [], "B": []}
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or line[76:78].strip() in {"H", "D"}:
            continue
        chain = line[21]
        if chain in atoms:
            xyz = [float(line[start : start + 8]) for start in (30, 38, 46)]
            atoms[chain].append((line[17:20], line[22:27], xyz))
    return atoms


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    args = parser.parse_args()

    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(
        ["complex", "ligand_heavy_atoms", "contact_ligand_atoms_4.5A", "contact_residues_4.5A",
         "min_distance_A", "pairs_under_2A", "cys_residues_within_4.5A"]
    )
    found = 0
    for path in sorted(args.run_dir.glob("job_*/*.pdb")):
        if path.stem.endswith("_binder"):
            continue
        atoms = read_heavy_atoms(path)
        ligand, protein = atoms["A"], atoms["B"]
        if not ligand or not protein:
            raise SystemExit(f"Missing ligand chain A or binder chain B in {path}")
        lig_xyz = np.array([atom[2] for atom in ligand])
        prot_xyz = np.array([atom[2] for atom in protein])
        distances = np.linalg.norm(lig_xyz[:, None] - prot_xyz, axis=2)
        near = distances <= 4.5
        contact_residues = {protein[i][:2] for i in np.flatnonzero(near.any(axis=0))}
        writer.writerow([
            path.stem, len(ligand), int(near.any(axis=1).sum()), len(contact_residues),
            f"{distances.min():.2f}", int((distances < 2.0).sum()),
            sum(name == "CYS" for name, _ in contact_residues),
        ])
        found += 1
    if not found:
        raise SystemExit(f"No generated complexes found in {args.run_dir}")


if __name__ == "__main__":
    main()
