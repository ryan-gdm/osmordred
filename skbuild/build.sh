#!/bin/bash
set -e

export RDKIT_VERSION=2023.9.3
export PYTHON_VERSION=3.11
export BOOST_VERSION=1.82.0

echo "Removing existing environment (if present)"
conda env remove -y -n osmordred || true

echo "Creating new conda environment"
conda create -y -n osmordred \
    python=${PYTHON_VERSION} \
    boost=${BOOST_VERSION} \
    eigen \
    lapack \
    ninja \
    python-build \
    rdkit=${RDKIT_VERSION} \
    blas=*=*mkl \
    -c conda-forge

echo "Cleaning previous build artifacts"
rm -rf _skbuild osmordred.egg-info dist

# Initialize conda for bash shell
eval "$(conda shell.bash hook)"

echo "Activating the environment"
conda activate osmordred

echo "Building the wheel"
conda run -n osmordred python -m build

echo "✅ Wheel build complete"

echo "Installing the wheel"
pip install dist/osmordred-0.2.0-cp311-cp311-linux_x86_64.whl --force-reinstall

echo "✅ Wheel installation complete"
