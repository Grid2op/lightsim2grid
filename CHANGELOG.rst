Change Log
===========
[TODO]
--------
- improve speed
- improve documentation
- configure circleci
- configure read the docs
- easier building (get rid of the "make" part)
- build on windows
- distribute pypi package
- try to implement an easier powerflow (for example gauss siedel) to make it available on all platform
- make dc approx available on all platform

[0.3.0] - 2020-10-xx
-------------------------
- [ADDING] Support for pickle for the lightsim Backend.
- [ADDING] LightSim should now be compatible with windows (implementation of a powerflow mode without
  using the SuiteSparse KLU linear solver but rather the Eigen SparseLU one)

[0.2.4] - 2020-08-xx
--------------------
- [FIXED] issue for copying environment

[0.2.3] - 2020-08-03
--------------------
- [UPDATED] consistent behaviour between grid2op.PandaPowerBackend and LightSimBackend for action that
  set the bus of only one extremity of a powerline.
- [ADDED] compatibility with grid2op 1.2.0

[0.2.2] - 2020-06-25
---------------------
- [UPDATED] removing the `-march=native` that causes some difficulty for some compilers
- [ADDED] compatibility with grid2op 1.0.0

[0.2.1] - 2020-06-xx
--------------------
- [FIXED] update of the `topo_vect` attribute in class `LightSimBackend` when reset.
- [ADDED] a github issue template

[0.2.0] - 2020-06-15
--------------------
- [ADDED] the changelog
- [FIXED] the import of files when elements where not in service
- [FIXED] a bad catch of a divergence in the solver
- [IMPROVED] the speed to apply the actions
- [FIXED] tests for the backend in grid2op and here are not identical without (too much) duplicates
