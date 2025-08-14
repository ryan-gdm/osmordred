#!/bin/bash
set -e

# === CONFIG ===
export PYTHON_VERSION=3.12
export BOOST_VERSION="1.86.*"
# ==============

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
    rdkit-dev \
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

# Detect the correct cpXXX tag for the installed Python
PYTHON_TAG=$(python -c "import sys; print(f'cp{sys.version_info.major}{sys.version_info.minor}')")

# Find matching wheel in dist/
WHEEL_FILE=$(ls dist/*${PYTHON_TAG}*.whl | head -n 1)
if [[ -z "$WHEEL_FILE" ]]; then
    echo "❌ Could not find wheel for ${PYTHON_TAG}"
    exit 1
fi

echo "Installing the wheel: $WHEEL_FILE"
pip install "$WHEEL_FILE" --force-reinstall

echo "✅ Wheel installation complete"
echo "Running post-installation tests..."

cd test
python test.py && echo "✅ Tests passed" || { echo "❌ Tests failed"; exit 1; }
