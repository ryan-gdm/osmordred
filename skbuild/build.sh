#!/bin/bash
set -e

# clean up generated dirs
rm -rf _skbuild src/cppmordred.egg-info dist

python setup.py build
python setup.py sdist bdist_wheel
echo Wheel build complete