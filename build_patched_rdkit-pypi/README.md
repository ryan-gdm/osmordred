# Building patched rdkit wheel
Run the script with the necessary CI Buildwheel env vars, e.g.:
```
CIBW_PLATFORM=linux CIBW_BUILD=cp311-manylinux_x86_64 ./build_rdkit-pypi.sh
```