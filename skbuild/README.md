python installation:

./build.sh

conda activate cpp_mordred

pip install dist/cppmordred-0.1.0-cp311-cp311-macosx_15_0_arm64.whl --force-reinstall

# note that attempting to `import cppmordred` from python running in the skbuild directory does not work, as it imports the `cppmordred` subdirectory.
cd test
python test.py
