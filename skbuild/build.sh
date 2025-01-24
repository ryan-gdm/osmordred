#!/bin/bash
set -e

# clean up generated dirs
rm -rf _skbuild src/cppmordred.egg-info dist

conda run -n cpp_mordred python -m build

echo Wheel build complete