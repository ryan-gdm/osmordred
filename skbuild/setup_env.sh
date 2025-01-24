#!/bin/bash
set -e

eval "$(conda shell.bash hook)"

conda create -y -n "cpp_mordred" python=3.11
conda activate cpp_mordred

conda install -y boost==1.82.0 eigen lapack mkl ninja rdkit==2023.9.3 -c conda-forge

pip install scikit-build tqdm
