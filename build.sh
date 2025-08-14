#!/bin/bash
set -e

export RDKIT_VERSION=2025.03.1
export PYTHON_VERSION=3.13
export BOOST_VERSION="1.86.*"

echo "Removing existing environment (if present)"
conda env remove -y -n osmordred || true

echo "Creating new conda environment"
conda create -y -n osmordred \
    python=${PYTHON_VERSION} \
    libboost-python=${BOOST_VERSION} \
    libboost-devel=${BOOST_VERSION} \
    libboost=${BOOST_VERSION} \
    eigen \
    lapack \
    ninja \
    python-build \
    rdkit-dev=${RDKIT_VERSION} \
    blas=*=*mkl \
    -c conda-forge

echo "Cleaning previous build artifacts"
rm -rf _src osmordred.egg-info dist

# Initialize conda for bash shell
eval "$(conda shell.bash hook)"

echo "Activating the environment"
conda activate osmordred

echo "Building the wheel"
export CMAKE_PREFIX_PATH=${CONDA_PREFIX}
conda run -n osmordred python -m build

echo "✅ Wheel build complete"

echo "Installing the wheel"
pip install dist/osmordred-0.3.0-cp313-cp313-linux_x86_64.whl --force-reinstall

echo "✅ Wheel installation complete"
echo "Running post-installation tests..."

cd test
python test.py && echo "✅ Tests passed" || { echo "❌ Tests failed"; exit 1; }