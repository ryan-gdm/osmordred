#!/bin/bash
set -e

function print_error {
    echo
    echo "Failed to create conda environment!!"
    exit 1
}
trap print_error ERR

echo "Removing existing environment (if present)"
conda env remove -y -n cpp_mordred &>/dev/null || true

conda_packages="boost==1.82.0 eigen lapack ninja python-build rdkit==2023.9.3"
if [[ "$OSTYPE" =~ ^darwin.* ]]; then
    echo "Creating conda env with MacOS packages"
    conda_packages="$conda_packages blas=*=*openblas"
elif [[ "$OSTYPE" =~ ^linux.* ]]; then
    echo "Creating conda env with Linux packages"
    conda_packages="$conda_packages blas=*=*mkl"
else
    echo "Don't recogize os: $OSTYPE"
    exit 1
fi

conda create -y -n cpp_mordred $conda_packages python=3.11 -c conda-forge
