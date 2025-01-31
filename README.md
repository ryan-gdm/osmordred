This repo was inspired by the Dragon, Padel, Mordred and other toolkits to generate empirical molecular features

The intention is to get 2D major molecular descriptors fused with rdkit backend at cpp level to get very fast computation in parallel if needed.


Remark on reproductibility:

The author spent quite time to implement a descent Information Content descriptor version based on the first paper of Basak in 1984 idea.
Information Content is not 100% identical but close enought to the logic to mimic the Information Content features.
This was indeed during this period that Triplet features from Basak team was also implemented. 

Future:
Current version is around 10k lines of codes in only one file. It will be refactor to follow standard RDkit python bindings.
Additionally a list of other descriptors were added to produce now 3586 individual features.
Code can be run on linux and mac (tested with ARM mac and Ubuntu 22.04 version) using python scripts or functions call.


Speed:
This is fully parallelized. Lapack was selected of it speed specially on the SVD decompostion of symmetrical squared matrix instead of Eigen3 solvers. 
The Lapack can produce very small fluctuation for almost zeros Eigen values and affect very slighlty few descriptors.   


Installation:

Method 1 : from scratch create a new environment python 3.11

./setup_env.sh

./build.sh

conda activate cpp_mordred

pip install dist/cppmordred-0.1.0-cp311-cp311-macosx_15_0_arm64.whl --force-reinstall or linux wheel file


Method 2 : include into your current environement (for the moment python 3.11 for RDKit 2023.9.3)

rm -rf _skbuild src/cppmordred.egg-info dist
- on mac: conda install boost==1.82.0 eigen lapack ninja python-build rdkit==2023.9.3 blas='*=*openblas' -c conda-forge
- on linux: conda install boost==1.82.0 eigen lapack ninja python-build rdkit==2023.9.3 blas='*=*mkl' -c conda-forge
conda run python -m build

Final installation of the wheel (linux or mac):
pip install dist/cppmordred-0.1.0-cp311-cp311-linux_x86_64.whl --force-reinstall 
or 
pip install dist/cppmordred-0.1.0-cp311-cp311-macosx_15_0_arm64.whl --force-reinstall

# normally you can see the installation in your environement 
pip show cppmordred  

# note that attempting to `import cppmordred` from python running in the skbuild directory does not work, as it imports the `cppmordred` subdirectory.

Testing:
cd test
python test.py
