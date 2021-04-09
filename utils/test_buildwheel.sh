#!/bin/bash
# what i did manually
# docker run -t -d --name manylinux_wheel -v $PWD:/lightsim2grid quay.io/pypa/manylinux1_x86_64
# docker exec -it manylinux_wheel bash
# cd /lightsim2grid
# make
# cd ..
# /opt/python/cp37-cp37m/bin/pip wheel /lightsim2grid/ --no-deps -w wheelhouse/
# auditwheel repair wheelhouse/LightSim2Grid-0.4.0-cp37-cp37m-linux_x86_64.whl --plat manylinux1_x86_64 -w /lightsim2grid/dist/

# skip "manylinux1_i686" because of a "PRE_CMD"
# cd ..  # TODO change that to be more consistant and less system dependant
for PLAT in manylinux1_x86_64 manylinux2010_x86_64; do
    # update docker image
    DOCKER_IMAGE=quay.io/pypa/$PLAT
    docker pull $DOCKER_IMAGE
    # some "logs"
    echo "building wheel for "$PLAT
    echo "Looking for lightsim2grid source file in "`pwd`
    # run the command to build the wheels
    docker run --rm -e PLAT=$PLAT -v `pwd`:/lightsim2grid $DOCKER_IMAGE /lightsim2grid/utils/build-wheels.sh
done
