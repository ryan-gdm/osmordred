#!/bin/bash
set -e

# clean up generated dirs
rm -rf _skbuild src/cppmordred.egg-info dist

conda run -n cpp_mordred python setup.py build
conda run -n cpp_mordred python setup.py sdist bdist_wheel
echo Wheel build complete