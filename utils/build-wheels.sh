#!/bin/bash
set -e -u -x

# largely inspired from https://github.com/pypa/python-manylinux-demo/blob/master/travis/build-wheels.sh

function repair_wheel {
    wheel="$1"
    if ! auditwheel show "$wheel"; then
        echo "Skipping non-platform wheel $wheel"
    else
        auditwheel repair "$wheel" --plat "$PLAT" -w /lightsim2grid/wheelhouse/
    fi
}

# see https://github.com/pypa/manylinux/issues/64#issuecomment-220997029
# for restraining python version
PYBINS=(
  "/opt/python/cp36-cp36m/bin"
  "/opt/python/cp37-cp37m/bin"
  "/opt/python/cp38-cp38m/bin"
  "/opt/python/cp39-cp39m/bin"
  )
# see https://github.com/pypa/manylinux/issues/64#issuecomment-221022231 for another alternative

# create the static files for SuiteSparse
cd /lightsim2grid
make clean
make
cd ..

# create the wheelhouse file where the... wheels will be housed XD
mkdir wheelhouse

# Compile wheels of the python package
#for PYBIN in ${PYBINS[@]}; do
for PYBIN in /opt/python/cp3*/bin; do
    "${PYBIN}/pip" install pybind11
    "${PYBIN}/pip" wheel /lightsim2grid/ --no-deps -w wheelhouse/
done

# clean the suisparse repo after making wheels
cd /lightsim2grid
make clean
cd ..


# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    repair_wheel "$whl"
done

# Install packages and test
# TODO
# for PYBIN in ${PYBINS[@]}; do
#    "${PYBIN}/pip" install lightsim2grid --no-index -f /lightsim2grid/wheelhouse
    # (cd "$HOME"; "${PYBIN}/nosetests" pymanylinuxdemo)  # TODO add the tests
# done
