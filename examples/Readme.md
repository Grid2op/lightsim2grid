# Code examples

In this folder are provided some example code using lightsim2grid.

Feel free to reuse these examples to match your needs.

## Description

- `contingency_analysis` shows how to use the `ContingencyAnalysis` module to simulate disconnections of powerline one after the other and benchmarks grid2op raw usage. It also check that that the results are the same when compared with grid2op.
- `time_serie` shows how to use the `TimeSerie` module efficiently and benchmarks this module compared to raw grid2op usage. It also check that that the results are the same when compared with grid2op.
- (advanced) `timeseries_with_grid2op` details how to use the c++ `TimeSeries` class
- (advanced) `timeseries_with_grid2op_multithreading` explains how to use the c++ `TimeSeries` class 
  in parrallel settings
- (advanced) `external_algorithm` a minimal, from-scratch solver plugin (`DummyExternalAlgo`), showing
  the whole plugin mechanism: building against an installed lightsim2grid with CMake,
  `load_algorithm_plugin(...)`, and selecting it with `change_algorithm("DummyExternal")`.
  See the "External Algorithm Plugins" page of the documentation for the full write-up.
- (advanced) `dist_slack_algorithm` a plugin registering a distributed-slack Newton-Raphson variant
  (`NRAlgoDistSlack`) with both a `SparseLU` and a `KLU` linear-solver backend.
- (advanced) `lm_algorithm` a plugin registering a Levenberg-Marquardt-damped Newton-Raphson variant
  (`LMNRAlgo`), again with `SparseLU` and `KLU` backends.
- `cmake/MatchLightsim2gridBuildFlags.cmake` a CMake helper, included by the three plugin examples
  above, that matches a plugin's `-march=native` / `-O3` flags to the ones `lightsim2grid_core` was
  built with -- necessary to avoid Eigen alignment mismatches between the two shared libraries.
