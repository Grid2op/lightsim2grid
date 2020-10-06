Getting started
===================================

In this chapter we present how to install lightsim2grid.

Installation
############

.. _requirements:

Requirements
*************
The requirements are:

- a working compiler
- python >= 3.6
- git
- the python packages "pybind11"

Install a compiler
~~~~~~~~~~~~~~~~~~~
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
feel free to send us a `github issue <https://github.com/BDonnot/lightsim2grid/issues>`_ and we will
do our best to answer.

Install python and git
~~~~~~~~~~~~~~~~~~~~~~~~~
Once you have a compiler you need to install **python** (again we will not cover how to get python on your
system) [because this is a python package] and **git** to install this package easily.

.. _install pybind11:

Install "pybind11"
~~~~~~~~~~~~~~~~~~~~~~~~
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


Installation
*************

This package is not (yet) available on pypi. Which means you have to install it from sources.

Do not forget to look at the :ref:`requirements` before proceeding in this step.

1. Retrieve the sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, you can download it with git with:

.. code-block:: bash

    git clone https://github.com/BDonnot/lightsim2grid.git
    cd lightsim2grid
    # it is recommended to do a python virtual environment
    python -m virtualenv venv  # optional
    source venv/bin/activate  # optional

    # retrieve the code of SparseSuite and Eigen (dependencies, mandatory)
    git submodule init
    git submodule update


Compilation of SuiteSparse (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SuiteSparse comes with the faster KLU linear solver.

Since version 0.3.0 this requirement has been removed. This entails
that on linux / macos you can still benefit from the faster KLU solver. On windows you will still benefit from the
speed up of lightsim (versus the default PandaPowerBackend) but this speed up will be less than if you manage
to compile SuiteSparse (see the subsection [Benchmark](#benchmark) for more information).

**NB** in both cases the algorithm to compute the powerflow is exactly the same. It is a
Newton Raphson based method. But to carry out this algorithm, one need to solver some linear equations. The only
difference in the two version (with KLU and without) is that the linear equation solver is different. Up to the
double float precision, both results (with and without KLU) should match.

We only detail the compilation on a system using "make" (so most likely GNU-Linux and MacOS). If you manage to
do this step on Windows, you can continue (and let us know!). If you don't feel confortable with this, we
provided a docker version. See the next section for more information.


.. code-block:: bash

    # compile static libraries of SparseSuite
    make

And yes that is it :-)

2. Installation of the python package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Now you simply need to install the lightsim2grid package this way, like any python package:

.. code-block:: bash

    # compile and install the python package
    python3 -m pip install -U .

**NB** please refer to the section :ref:`install pybind11` for more information. Indeed the
command to invoke python may vary. You may need to replace `python3` with `python`, `py` or `py3` for
example.

And you are done :-)

Start Using LightSim2grid
############################
The preferred way to use light2im simulator is with Grid2op. And in this case, you can simply use it
this way:

.. code-block:: python

    import grid2op
    from lightsim2grid.LightSimBackend import LightSimBackend
    from grid2op.Agent import RandomAgent

    # create an environment
    env_name = "rte_case14_realistic"  # for example, other environments might be usable
    env = grid2op.make(env_name,
                       backend=LightSimBackend()  # this is the only change you have to make!
                       )

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


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`