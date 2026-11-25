#!/usr/bin/env bash
set -eo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
work=$(dirname "$repo")
cd "$work"
exec > >(tee "$work/pilot.log") 2>&1
export XDG_CACHE_HOME="$work/cache"
export XDG_CONFIG_HOME="$work/config"
export XDG_DATA_HOME="$work/data"
export UV_PYTHON_INSTALL_DIR="$work/python"
export UV_INSTALL_DIR="$work/bin"
export PATH="$UV_INSTALL_DIR:$PATH"
if ! command -v uv >/dev/null; then
    curl -fsSL https://astral.sh/uv/install.sh | env UV_NO_MODIFY_PATH=1 sh
fi
if command -v nvidia-smi >/dev/null; then
    nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv
fi
if [[ ! -d Proteina-Complexa/.git ]]; then
    git clone --depth 1 https://github.com/NVIDIA-BioNeMo/Proteina-Complexa.git Proteina-Complexa
fi
cd Proteina-Complexa
revision=618a5a0b34a5b038e8c124315f395273d71537ff
if ! git cat-file -e "$revision^{commit}"; then
    git fetch --depth 1 origin "$revision"
fi
git checkout --detach "$revision"
if [[ ! -f .pilot-install-complete ]]; then
    bash env/build_uv_env.sh --minimal
    touch .pilot-install-complete
fi
source .venv/bin/activate
python -c 'import torch; print("CUDA:", torch.cuda.get_device_name())'
if [[ ! -f .env ]]; then
    complexa init
    sed -i "s|^LOCAL_CODE_PATH=.*|LOCAL_CODE_PATH=$PWD|; s|^LOCAL_DATA_PATH=.*|LOCAL_DATA_PATH=$PWD/assets|" .env
fi
if [[ ! -f env.sh ]]; then
    complexa init uv
fi
source env.sh
mkdir -p ckpts
weights_url='https://api.ngc.nvidia.com/v2/models/org/nvidia/team/clara/proteina_complexa_ligand/1.0/files?redirect=true&path='
for checkpoint in complexa_ligand.ckpt complexa_ligand_ae.ckpt; do
    if python -c 'import sys, zipfile; sys.exit(not zipfile.is_zipfile(sys.argv[1]))' "ckpts/$checkpoint"; then
        continue
    fi
    temporary=$(mktemp "ckpts/$checkpoint.XXXXXX")
    curl --fail --location --retry 3 --silent --show-error --output "$temporary" "$weights_url$checkpoint"
    python -c 'import sys, zipfile; sys.exit(not zipfile.is_zipfile(sys.argv[1]))' "$temporary"
    mv "$temporary" "ckpts/$checkpoint"
done
bash "$repo/models/proteina-complexa/run_pilot.sh" "$PWD"
