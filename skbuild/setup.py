from skbuild import setup
from setuptools import find_packages
import sysconfig
import os



setup(
    name="cppmordred",
    version="0.1.3",
    description="A Python package for mordred in cpp using RDKit 2023.09.3, Lapack.",

    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Guillaume Godin",
    author_email="guillaume.godin@gmail.com",

    license="BSD-3-Clause",
    packages=["cppmordred"],
    cmake_args=[
        f"-DCMAKE_PREFIX_PATH={os.environ.get('CONDA_PREFIX')}",
    ],
    python_requires=">=3.11.8",
    install_requires=["scikit-build", "numpy==1.26.4"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD-3-Clause License",
        "Operating System :: OS Independent",
    ],
)
