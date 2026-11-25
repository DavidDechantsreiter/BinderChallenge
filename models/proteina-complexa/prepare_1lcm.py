"""Prepare a heavy-atom copy of 1LCM with the missing FGA OXT restored."""

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
from biotite.structure.io.pdbx import CIFCategory, CIFFile


def unit(vector):
    return vector / np.linalg.norm(vector)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    cif = CIFFile.read(args.source)
    sites = cif.block["atom_site"]
    columns = {name: sites[name].as_array(str) for name in sites}
    heavy = ~np.isin(columns["type_symbol"], ["H", "D"])
    if Counter(columns["type_symbol"][heavy]) != {"C": 49, "N": 10, "O": 11}:
        raise ValueError("Expected the original 1LCM heavy-atom composition.")

    glutamate = (columns["auth_asym_id"] == "A") & (columns["auth_seq_id"] == "6")
    xyz_fields = ("Cartn_x", "Cartn_y", "Cartn_z")
    coordinates = np.column_stack([columns[field].astype(float) for field in xyz_fields])

    def atom_index(name):
        indices = np.flatnonzero(glutamate & (columns["auth_atom_id"] == name))
        if len(indices) != 1:
            raise ValueError(f"Expected exactly one FGA atom named {name}.")
        return indices[0]

    carbon = atom_index("C")
    oxygen = atom_index("O")
    alpha = atom_index("CA")
    center = coordinates[carbon]
    direction = -unit(unit(coordinates[alpha] - center) + unit(coordinates[oxygen] - center))
    added = center + 1.30 * direction

    # Insert into the existing residue. Preserve the CIF chemical component
    # dictionary and cross-residue bonds so the loader can recover bond orders.
    template = {name: values[oxygen] for name, values in columns.items()}
    template.update(auth_atom_id="OXT", label_atom_id="OXT")
    template.update({field: f"{value:.3f}" for field, value in zip(xyz_fields, added)})
    insertion = int(np.count_nonzero(heavy[: oxygen + 1]))
    revised = {
        name: np.insert(values[heavy].astype(object), insertion, template[name]).tolist()
        for name, values in columns.items()
    }
    revised["id"] = [str(i) for i in range(1, len(revised["id"]) + 1)]
    cif.block["atom_site"] = CIFCategory(revised)
    cif.write(args.output)

    heavy[[carbon, oxygen, alpha]] = False
    closest = min(np.linalg.norm(coordinates[heavy] - added, axis=1))
    print(f"Added FGA OXT (C–OXT = 1.30 Å, nearest other heavy atom {closest:.2f} Å).")


if __name__ == "__main__":
    main()
