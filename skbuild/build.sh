#!/bin/bash
set -e

# clean up generated dirs
rm -rf _skbuild src/osmordred.egg-info dist

conda run -n osmordred python -m build

echo Wheel build complete
