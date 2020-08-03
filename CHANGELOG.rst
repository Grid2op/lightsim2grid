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
