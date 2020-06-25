# Copyright (c) 2019-2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import os
import argparse
import re
import pdb
import subprocess
import time


def start_subprocess_print(li, sleepbefore=2, cwd=None):
    print("Will execute command after {}s: \n\t{}".format(sleepbefore, " ".join(li)))
    time.sleep(sleepbefore)
    subprocess.run(li, cwd=cwd)


def modify_and_push_docker(path,
                           ls_version,
                           docker_tags=[],
                           docker_extra_build_arg=[],
                           dockerfile_name="Dockerfile"):
    # Create new docker containers
    li_tags = []
    for vers_ in docker_tags:
        li_tags += ["-t", "{}/lightsim2grid:{}".format(dockeruser, vers_)]
    start_subprocess_print(["docker", "build"] + docker_extra_build_arg +
                           li_tags +
                           ["--build-arg", "ls_version=v{}".format(ls_version)] +
                           ["-f", "{}".format(os.path.join(path, dockerfile_name)), "."],
                           cwd=path)

    # push the containers
    for vers_ in docker_tags:
        start_subprocess_print(["docker", "push", "{}/lightsim2grid:{}".format(dockeruser, vers_)], cwd=path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Update the version of grid2op in the python files.')
    parser.add_argument('--version', default=None,
                        help='The new version to update.')
    parser.add_argument('--dockeruser', default='bdonnot',
                        help='The docker hub username.')
    parser.add_argument('--path', default=os.path.abspath("./utils"),
                        help='The path of the root directory of Grid2op (default {}'.format(os.path.abspath(".")))
    args = parser.parse_args()
    path = args.path
    dockeruser = args.dockeruser
    version = args.version

    if args.version is None:
        raise RuntimeError("script \"update_version\" should be called with a version number.")

    try:
        maj_, min_, minmin_, *post = version.split(".")
    except:
        raise RuntimeError(
            "script \"push_docker\": version should be formated as XX.YY.ZZ (eg 0.3.1). Please modify \"--version\" "
            "argument")

    regex_version = "[0-9]+\.[0-9]+\.[0-9]+(.post[0-9]+){0,1}"
    if re.match("^{}$".format(regex_version), version) is None:
        raise RuntimeError(
            "script \"push_docker\": version should be formated as XX.YY.ZZ (eg 0.3.1) and not {}. Please modify "
            "\"--version\" argument".format(version))

    modify_and_push_docker(path=path,
                           ls_version=version,
                           docker_tags=["test"],
                           dockerfile_name="Dockerfile_test",
                           docker_extra_build_arg=["--no-cache"])

    modify_and_push_docker(path=path,
                           ls_version=version,
                           docker_tags=["light"],
                           dockerfile_name="Dockerfile_light",
                           docker_extra_build_arg=["--no-cache"])

    modify_and_push_docker(path=path,
                           ls_version=version,
                           docker_tags=["{}".format(version), "latest"],
                           dockerfile_name="Dockerfile")
