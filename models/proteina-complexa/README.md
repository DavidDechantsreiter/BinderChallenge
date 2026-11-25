# Proteina-Complexa pilot with microcystin-LR

Proteina-Complexa generated two complexes for 1LCM, each with a 100-residue protein and all 71 ligand heavy atoms. Both predicted poses have five interatomic pairs below 2 Å, including contacts with backbone oxygens.

| Design | Ligand atoms within 4.5 Å | Binder residues within 4.5 Å | Closest pair | Pairs below 2 Å | Cysteines within 4.5 Å |
|---|---:|---:|---:|---:|---:|
| `job_0_n_100_id_0_single_orig0` | 50 of 71 | 19 | 1.36 Å | 5 | 0 |
| `job_0_n_100_id_1_single_orig0` | 53 of 71 | 22 | 1.50 Å | 5 | 0 |

See the [structural analysis and animations](ANALYSIS.md) for the original PDBs, clash locations, tertiary folds, and approximate volumes.

## Input and execution

The original `1LCM.pdb` contains 70 heavy atoms. It lacks the `OXT` oxygen of glutamate `FGA` 6. `prepare_1lcm.py` added that oxygen to a copy of the mmCIF at an estimated position 1.30 Å from the carbon. The 70 existing heavy atoms retain the original PDB coordinates. The Complexa loader verified that all 71 heavy atoms were connected and that the ring closure and bond orders were present.

Job `proteina-recheck` completed on **25 November 2026** on a **40 GB A100** in the Run:ai project `make-dlai-romero`. The generation step took **75.4 s**. The pilot used NVIDIA revision `618a5a0b34a5b038e8c124315f395273d71537ff`, the ligand checkpoint, and two single-pass samples with batch size 1. The setup script checks the checkpoint archive before reuse and downloads replacements into a temporary file. [NVIDIA documentation](https://github.com/NVIDIA-BioNeMo/Proteina-Complexa/tree/618a5a0b34a5b038e8c124315f395273d71537ff)

The RCP clone is at `/rcp-home/binder-complexa-pilot/BinderChallenge`. The run directory contains the original PDBs and `pilot_geometry.csv`.

```text
/rcp-home/binder-complexa-pilot/Proteina-Complexa/inference/1lcm_pilot_20261125T073924Z/
```

To rerun from the clone in a GPU job with your personal volume mounted at `/rcp-home`, run:

```bash
bash models/proteina-complexa/setup_rcp.sh
```
