#!/bin/bash
set -e

# work in the script directory
cd $(dirname $0)

if [ -z "$CIBW_BUILD" ]; then
    echo "You must set CIBW_BUILD to set a build target environment."
    exit 1
fi

CONDA_ENV=build-rdkit-pypi
CONDA_RUN="conda run -n $CONDA_ENV"

echo "Will create a conda env ($CONDA_ENV) to launch cibuildwheel from"
conda create -y -n $CONDA_ENV python=3.11 
$CONDA_RUN pip install cibuildwheel

# there's no tag in rdkit-pypi for the 2023.9.3 release
# (I think it moved from a different repo)
RDKIT_PYPI_2023_9_3_SHA=311157810e3d018bd2333f0f3c75bcc8538bf486
git clone https://github.com/kuelumbus/rdkit-pypi.git
cd rdkit-pypi
git checkout -b 2023.9.3-osmordred $RDKIT_PYPI_2023_9_3_SHA
# Patch the setup.py and pyproject.toml
git apply ../rdkit-pypi_2023_09_3_osmordred.diff

cp -a ../../osmordred_rdkit_2023_09_3_patches .
mkdir osmordred_source
cp -a ../../Code osmordred_source

$CONDA_RUN python3 -m cibuildwheel --platform linux --output-dir wheelhouse --config-file pyproject.toml