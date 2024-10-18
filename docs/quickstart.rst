Getting started
===================================

In this chapter we present how to install lightsim2grid.

.. versionadded:: 0.5.4

    lightsim2grid can be installed directly from pypi.

Standard recommended installation
------------------------------------

Simply install lightsim2grid as any python package:

.. code-block:: bash
    
    pip install lightsim2grid

.. note::
    depending on your settings, OS, installation of python, you might want to do `python -m pip XXX`, 
    `python3 -m pip install XXX` or `py -m pip install XXX` etc. see the documentation of `pip` for your system)


Start Using LightSim2grid
---------------------------

The preferred way to use light2im simulator is with Grid2op. And in this case, you can simply use it
this way:

.. code-block:: python

    import grid2op
    from lightsim2grid.LightSimBackend import LightSimBackend
    from grid2op.Agent import RandomAgent

    # create an environment
    env_name = "l2rpn_case14_sandbox"  # for example, other environments might be usable
    env = grid2op.make(env_name,
                       backend=LightSimBackend()  # this is the only change you have to make!
                       )

    # As of now, you can do whatever you want with this grid2op environment
    # for example...

    # create an agent
    my_agent = RandomAgent(env.action_space)

    # proceed as you would any open ai gym loop
    nb_episode = 10
    for _ in range(nb_episde):
        # you perform in this case 10 different episodes
        obs = env.reset()
        reward = env.reward_range[0]
        done = False
        while not done:
            # here you loop on the time steps: at each step your agent receive an observation
            # takes an action
            # and the environment computes the next observation that will be used at the next step.
            act = agent.act(obs, reward, done)
            obs, reward, done, info = env.step(act)
            # the `LightSimBackend` will be used to carry out the powerflow computation instead
            # of the default grid2op `PandaPowerBackend`


Installation from source (advanced usage)
-------------------------------------------

To install this package from source, you need, in summary, to:

- clone this repository and get the code of Eigen (mandatory for compilation) and SparseSuite (optional, but recommended)
- (optional, but recommended) compile a piece of SparseSuite
- (optional) [experimental] retrieve and get a proper license for the NICSLU linear solver (see https://github.com/chenxm1986/nicslu)
- (optional) specify some compilation flags to make the package run faster on your machine
- install the package
  

.. _requirements:

Requirements
~~~~~~~~~~~~~~~~~~~~~~

The requirements are:

- a working compiler
- python >= 3.6
- git
- the python packages "pybind11"

Install a compiler
+++++++++++++++++++
It relies on c++ to carry out some computations faster than pure python solvers. To integrate
c++ into python the excellent `pybind11 <https://github.com/pybind/pybind11>`_ library is used.
This entails that you need to have a *compiler* that can be used by pybind11 (don't hesitate to
`check the list of supported compilers <https://github.com/pybind/pybind11#supported-compilers>`_
which was, at time of writing:

- Windows: Microsoft Visual Studio 2015 Update 3 or newer
  (see `here <https://visualstudio.microsoft.com/vs/features/cplusplus/>`_ for help on how to install it)
- Linux (Ubuntu, Fedora, etc.): GCC 4.8 or newer (on ubuntu `sudo apt install build-essential` or
  on Fedora: somthing like `sudo dnf install make automake gcc gcc-c++ kernel-devel`)
- MacOs: Clang/LLVM 3.3 or newer (for Apple Xcode’s clang, this is 5.0.0 or newer), you can
  install it by typing `brew install llvm` in a terminal.

We do not cover in this installation guide how to install such compiler. But if you have any issue,
feel free to send us a `github issue <https://github.com/Grid2Op/lightsim2grid/issues>`_ and we will
do our best to answer.

Install python and git
+++++++++++++++++++++++
Once you have a compiler you need to install **python** (again we will not cover how to get python on your
system) [because this is a python package] and **git** to install this package easily.

.. _install pybind11:

Install "pybind11"
++++++++++++++++++++++
As most python package you can install it using `pip`.

On MacOs or Linux (Ubuntu, Fedora, etc.) you can install it from a commandline with:

.. code-block:: bash

    python3 -m pip install pybind11

On windows, the commandline is harder to find, and the command to invoke python can vary (sometimes
it is `py` sometimes `python` etc.) depending on your installation. You can try:

.. code-block:: bash

    py -m pip install pybind11

or

.. code-block:: bash

    python3 -m pip install pybind11


Installation of the python package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that everything is setup, you can proceed with the installation of the package

1. Retrieve the sources
+++++++++++++++++++++++

First, you can download it with git with:

.. code-block:: bash

    git clone https://github.com/Grid2Op/lightsim2grid.git
    cd lightsim2grid
    # it is recommended to do a python virtual environment
    python -m virtualenv venv  # optional
    source venv/bin/activate  # optional

    # retrieve the code of SparseSuite and Eigen (dependencies, mandatory)
    git submodule init
    git submodule update


Compilation of SuiteSparse (optional)
+++++++++++++++++++++++++++++++++++++++

SuiteSparse comes with the faster KLU linear solver.

Since version 0.3.0 this requirement has been removed. This entails
that on linux / macos you can still benefit from the faster KLU solver. On windows you will still benefit from the
speed up of lightsim (versus the default PandaPowerBackend) but this speed up will be less than if you manage
to compile SuiteSparse (see the subsection *Benchmark* for more information).

**NB** in both cases the algorithm to compute the powerflow is exactly the same. It is a
Newton Raphson based method. But to carry out this algorithm, one need to solver some linear equations. The only
difference in the two version (with KLU and without) is that the linear equation solver is different. Up to the
double float precision, both results (with and without KLU) should match.

We only detail the compilation on a system using "make" (so most likely GNU-Linux and MacOS). If you manage to
do this step on Windows, you can continue (and let us know!). If you don't feel confortable with this, we
provided a docker version. See the next section for more information.

(optional) option A. Compilation of SuiteSparse using "make"
******************************************************************

This is the easiest method to compile SuiteSparse on your system but unfortunately it only works on OS where "make" is
available (*eg* Linux or MacOS) but this will not work on Windows... The compilation on windows is covered in the next
paragraph *(optional) option B. Compilation of SuiteSparse using "cmake"* bellow.

Anyway, in this case, it's super easy. Just do:

.. code-block:: bash

    # compile static libraries of SparseSuite
    make

And yes that is it :-)

(optional) option B. Compilation of SuiteSparse using "cmake"
******************************************************************
This works on most platform including MacOS, Linux and Windows.

It requires to install the free `cmake` program and to do a bit more works than for other system. This is why we
only recommend to use it on Windows.

The main steps (for windows, somme commands needs to be adapted on linux / macos) are:
1) `cd build_cmake`
2) `py generate_c_files.py`
3) `mkdir build` and cd there: `cd build`
4) `cmake -DCMAKE_INSTALL_PREFIX=..\built -DCMAKE_BUILD_TYPE=Release ..`
5) `cmake --build . --config Release`
6) `cmake --build . --config Release --target install`

For more information, feel free to read the dedicated [README](build_cmake/README.md).

(optional) Include NICSLU linear solver (experimental)
******************************************************************

Another linear solver that can be used with lighsim2grid is the "NICSLU" linear solver that might, in some cases, be
even faster than the KLU linear solver. This can lead to more speed up if using lighsim2grid.

To use it, you need to:

1) retrieve the sources (only available as a freeware) from https://github.com/chenxm1986/nicslu and save
   it on your machine. Say you clone this github repository in `NICSLU_GIT` 
   (*eg* NICSLU_GIT="/home/user/Documents/nicslu/"). Also note that you need to check that your usage
   is compliant with their license !
2) define the "PATH_NICSLU" environment variable **before** compiling lightsim2grid, on linux you can do
   `export PATH_NICSLU=NICSLU_GIT/nicsluDATE` 
   (for example `export PATH_NICSLU=/home/user/Documents/nicslu/nicslu202103` if you cloned the repository 
   as the example of `step 1)` and use the version of nicslu compiled by the author on March 2021 [version 
   distributed at time of writing the readme] )

And this is it. Lightsim will be able to use this linear solver.

Be carefull though, you require a license file in order to use it. As of now, the best way is to copy paste the license
file at the same location that the one you execute python from (*ie* you need to copy paste it each time).

.. _quickstart_install_custom:

(optional) customization of the installation
******************************************************************

If you bother to compile from source the package, you might also want to benefit from some
extra speed ups.

This can be achieve by specifying the `__O3_OPTIM` and `__COMPILE_MARCHNATIVE` environment variables. 

The first one will compile the package using the `-O3` compiler flag (`/O2` on windows) which will tell the compiler to optimize the code for speed even more.

The second one will compile the package using the `-march=native` flag (on macos and linux)

And example to do such things on a linux based machine is:

.. code-block:: bash

    export __O3_OPTIM=1
    export __COMPILE_MARCHNATIVE=1

On windows, the equivalent is either (*eg* in the good old fashion "cmd"):

.. code-block:: none

    setx __O3_OPTIM "1"
    setx __COMPILE_MARCHNATIVE "1"

or (*eg* in the "new" powershell):

.. code-block:: none

    $Env:__O3_OPTIM="1"
    $Env:__COMPILE_MARCHNATIVE="1"

If you want to disable them, you simply need to set their respective value to "0" instead of "1".

.. note::
    By default, on pypi, the packages are built with the `__O3_OPTIM` flag.

.. warning::
    If you use the `__COMPILE_MARCHNATIVE` flag, the python package you generate might not 
    work on other installation: you cannot reuse it on other computer, or on docker etc.

2. Installation of the python package
+++++++++++++++++++++++++++++++++++++++

Now you simply need to install the lightsim2grid package this way, like any python package:

.. code-block:: bash

    # compile and install the python package
    python3 -m pip install -U .

**NB** please refer to the section :ref:`install pybind11` for more information. Indeed the
command to invoke python may vary. You may need to replace `python3` with `python`, `py` or `py3` for
example.

.. note::
    We remind you that you can customize the installation, see more details in :ref:`quickstart_install_custom`

And you are done :-)

Usage with docker
-------------------------------------------

In this section we cover the use of docker with grid2op.

1. Install docker
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
First, you need to install docker. You can consult the 
`docker on windows <https://hub.docker.com/editions/community/docker-ce-desktop-windows>`_ if you use a windows like
operating system, if you are using MacOs you can consult 
`docker on Mac <https://hub.docker.com/editions/community/docker-ce-desktop-mac/>`_ . The installation of docker on linux
depends on your linux distribution, we will not list them all here.

2. Get the lightsim2grid image
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once done, you can simply "install" the lightsim2grid image with:

.. code-block:: bash

    docker pull bdonnot/lightsim2grid:latest


This step should be done only once (unless you delete the image) it will download approximately 4 or 5GB from the
internet. The lightsim2grid image contains lightsim and grid2op python packages (as well as their
dependencies), equivalent of what would be installed if you typed:
.. code-block:: bash

    pip install -U grid2op[optional] pybind11
    # and do steps detailed in section "Installation (from source)"
    # that we will not repeat


3. Run a code on this container
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can skip this section if you know how to use docker. We will present here "the simplest way" to use. This is NOT
a tutorial on docker, and you can find better use of this technology on 
`the docker website <https://www.docker.com/get-started>`_ .

For this tutorial, we suppose you have a script named `my_script.py` located in the directory (complete path) 
`DIR_PATH` (*e.g.* on windows you can have `DIR_PATH` looking like "c:\\User\\MyName\\L2RPNCompeitionCode" or 
on Linux `DIR_PATH` will look like "/home/MyName/L2RPNCompeitionCode", this path is your choice, you can name it
the way you like)

3.1) Start a docker container
+++++++++++++++++++++++++++++++++++++++

You first need to start a docker container and tell docker that the container can access your local files with the 
following command:

.. code-block:: bash

    docker run -t -d -p 8888:8888 --name lightsim_container -v DIR_PATH:/L2RPNCompeitionCode -w /L2RPNCompeitionCode bdonnot/lightsim2grid

More information on this command 
`in the official docker documentation <https://docs.docker.com/engine/reference/commandline/run/>`_

After this call you can check everything went smoothly with by invoking:

.. code-block:: bash

    docker ps

And the results should look like::

    CONTAINER ID        IMAGE                   COMMAND             CREATED             STATUS              PORTS               NAMES
    89750964ca55        bdonnot/lightsim2grid   "python3"           5 seconds ago       Up 4 seconds        80/tcp              lightsim_container


`DIR_PATH` should be replaced by the path on which you are working, see again the introduction of this
section for more information, in the example above this can look like:

.. code-block:: bash

    docker run -t -d -p 8888:8888 --name lightsim_container -v /home/MyName/L2RPNCompeitionCode:/L2RPNCompeitionCode -w /L2RPNCompeitionCode bdonnot/lightsim2grid


3.2) Execute your code on this container
+++++++++++++++++++++++++++++++++++++++++

Once everything is set-up you can execute anything you want on this container. Note that doing so, the execution
of the code will be totally independant of your system. Only the things located in `DIR_PATH` will be visible 
by your script, only the python package installed in the container will be usable, only the python interpreter
of the containter (python 3.6 at time of writing) will be usable etc.

.. code-block:: bash

    docker exec lightsim_container python my_script.py


Of course, the "my_script.py" should save its output somewhere on the hard drive.

If you rather want to execute a python REPL (read-eval-print loop), corresponding to the "interactive python 
interpreter", you can run this command:

.. code-block:: bash

    docker exec -it lightsim_container python


We also added the possibility to run jupyter notebook from this container. To do so, you can run the command:

.. code-block:: bash

    docker exec -it lightsim_container jupyter notebook --port=8888 --no-browser --ip='*' --allow-root

More information is provided in the official documentation of 
`docker exec  <https://docs.docker.com/engine/reference/commandline/exec/>`_.

3.3) Disclaimer
+++++++++++++++++++++++++++++++++++++++

Usually, docker run as root on your machine, be careful, you can do irreversible things with it. "A great power 
comes with a great responsibility".

Also, we recall that we presented a really short introduction to docker and its possibility. We have not implied
that this was enough, nor explain (on purpose, to make this short) any of the commands. 
We strongly encourage you to have a look for yourself. 

We want to recall the paragraph `7. Limitation of Liability` under which lightsim2grid, and this "tutorial" 
is distributed


.. note:: 
    
    Under no circumstances and under no legal 
    theory, whether tort (including negligence), 
    contract, or otherwise, shall any Contributor, or 
    anyone who distributes Covered Software as 
    permitted above, be liable to You for any direct, 
    indirect, special, incidental, or consequential 
    damages of any character including, without 
    limitation, damages for lost profits, loss of 
    goodwill, work stoppage, **computer failure** or
    **malfunction**, or any and all other commercial 
    damages or losses, even if such party shall have 
    been informed of the possibility of such damages.

1. Clean-up
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once you are done with your experiments, you can stop the docker container:

.. code-block:: bash

    docker container stop lightsim_container


This will free all the CPU / GPU resources that this container will use. If you want to start it again, for another 
experiment for example, just use the command:

.. code-block:: bash

    docker container start lightsim_container

This will allow you to run another batch of `dcoker exec` (see `3.2) Execute your code on this container`) 
without having to re run the container.


If you want to go a step further, you can also delete the container with the command:

.. code-block:: bash

    docker container rm lightsim_container

This will remove the container, and all your code executed there, the history of commands etc. If you want to use
lightsim2grid with docker again you will have to go through section `3. Run a code on this container` all over
again.

And if you also want to remove the image, you can do:

.. code-block:: bash

    docker rmi bdonnot/lightsim2grid 

**NB** this last command will completely erase the lightsim2grid image from your machine. This means that 
if you want to use it again, you will have to download it again (see section `2. Get the lightsim2grid image`)

Finally, you can see the official documentation in case you need to uninstall docker completely from your system.

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`