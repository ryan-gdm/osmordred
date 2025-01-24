from skbuild import setup
from setuptools import find_packages
import sysconfig
import os



setup(
    name="cppmordred",
    version="0.1.0",
    description="A Python package for mordred in cpp using RDKit 2023.09.3, LAPACK and pybind11.",

    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Guillaume Godin",
    author_email="guillaume.godin@gmail.com",

    license="MIT",
    packages=["cppmordred"],
    cmake_args=[
        # '-DCMAKE_MAKE_PROGRAM:STRING=make',
        # f"-DPython3_EXECUTABLE=/Users/guillaume-osmo/miniconda3/envs/osmo-sandox/bin/python",
        # f"-DRDKit_DIR=/Users/guillaume-osmo/miniconda3/pkgs/rdkit-2023.09.3-py311h64bc748_1/",
        # f"-DBOOST_ROOT=/Users/guillaume-osmo/miniconda3/envs/osmo-sandox/",
        f"-DLAPACK_LIBRARIES=/opt/homebrew/opt/lapack/lib",
        f"-DLAPACK_INCLUDE_DIRS=/opt/homebrew/opt/lapack/include",
        # f"-DPYTHON_SITE_PACKAGES=/Users/guillaume-osmo/miniconda3/envs/osmo-sandox/lib/python3.11/site-packages"
    ],
    python_requires=">=3.11.8",
    install_requires=["scikit-build", "numpy==1.26.4"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)