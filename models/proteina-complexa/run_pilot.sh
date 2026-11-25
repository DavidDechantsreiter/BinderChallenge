#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
source_cif="$repo/related_work/2026-microcystein/02_structures/microcystein/1LCM.cif"
complexa_dir=${1:?Usage: bash run_pilot.sh /path/to/Proteina-Complexa}
cd "$complexa_dir"
set +u
source .venv/bin/activate
source env.sh
set -u
run_name="1lcm_pilot_$(date -u +%Y%m%dT%H%M%SZ)"
input_dir="$PWD/inputs/$run_name"
mkdir -p "$PWD/inputs"
mkdir "$input_dir"
target="$input_dir/1LCM_completed.cif"
results="$PWD/inference/$run_name"
python "$repo/models/proteina-complexa/prepare_1lcm.py" "$source_cif" "$target"

# Check the molecular graph as seen by the loader used in Complexa.
python - "$target" <<'PY'
from collections import Counter
from sys import argv

import numpy as np
from atomworks.io.utils.io_utils import load_any

target = load_any(argv[1])[0]
heavy = target[~np.isin(np.char.upper(target.element), ["H", "D"])]
elements = Counter(np.char.upper(heavy.element))
expected = {"C": 49, "N": 10, "O": 12}
if elements != expected:
    raise SystemExit(f"Unexpected heavy-atom composition: {dict(elements)} != {expected}")

if heavy.bonds is None:
    raise SystemExit("The loader returned no heavy-atom bonds.")
if np.any(heavy.bonds.as_array()[:, 2] == 0):
    raise SystemExit("The loader returned untyped heavy-atom bonds.")
adjacency = heavy.bonds.adjacency_matrix()
seen = {0}
pending = [0]
while pending:
    for neighbor in np.flatnonzero(adjacency[pending.pop()]):
        if neighbor not in seen:
            seen.add(neighbor)
            pending.append(neighbor)
if len(seen) != len(heavy):
    raise SystemExit(f"Disconnected ligand: reached {len(seen)} of {len(heavy)} heavy atoms.")

first = np.flatnonzero((heavy.res_id == 1) & (heavy.atom_name == "N"))
last = np.flatnonzero((heavy.res_id == 7) & (heavy.atom_name == "C"))
if len(first) != 1 or len(last) != 1 or not adjacency[first[0], last[0]]:
    raise SystemExit("The loader did not preserve the bond closing the peptide ring.")
print(f"1LCM preflight passed: {len(heavy)} connected heavy atoms and intact ring closure.")
PY

# Hydra override selects the entire target rather than individual PDB residue names.
complexa generate --verbose configs/search_ligand_binder_local_pipeline.yaml \
    "++run_name=$run_name" \
    "++root_path='$results'" \
    ++generation.task_name=MC_LR_PILOT \
    "++generation.target_dict_cfg.MC_LR_PILOT={target_path:'$target',target_filename:'1LCM',source:local,pdb_id:'1lcm',binder_length:[100],hotspot_residues:[null],ligand:null,ligand_only:true,SMILES:null,use_bonds_from_file:true}" \
    ++generation.dataloader.batch_size=1 \
    ++generation.search.algorithm=single-pass \
    ++generation.reward_model=null \
    ++generation.dataloader.dataset.nres.nsamples=2

python "$repo/models/proteina-complexa/analyze_pilot.py" "$results" | tee "$results/pilot_geometry.csv"
