
Dependencies:
- rdkit 2023.9.3 from conda forge with include headers and lib
- eigen3 
- boost 1.82 (specially python)
- lapack (faster than eigen3)
- pybind11

# need to change pip rdkit to conda rdkit
- pip uninstall rdkit==2023.9.3 
- conda install rdkit==2023.9.3 -c conda-forge
- conda install boost==1.82.0 -c conda-forge

# adjust to your environment and also change the CMakelists.txt file accordingly
cmake -DPython3_EXECUTABLE=/Users/guillaume-osmo/miniconda3/envs/osmo/bin/python \
      -Dpybind11_DIR=$(python -m pybind11 --cmakedir) \
      -DRDKit_DIR=/Users/guillaume-osmo/miniconda3/pkgs/rdkit-2023.09.3-py311h64bc748_1/ \
      -DEigen3_DIR=/opt/homebrew/include/eigen3/  \
      ..




