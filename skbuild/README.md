python installation:

rm -rf _skbuild build

python setup.py build 

python setup.py sdist bdist_wheel

pip install dist/cppmordred-0.1.0-cp311-cp311-macosx_15_0_arm64.whl --force-reinstall  