Change Log
===========

[TODO]
--------
- [refacto] have a structure in cpp for the buses
- [refacto] have the id_grid_to_solver and id_solver_to_grid etc. directly in the solver and NOT in the gridmodel.
- [refacto] put some method in the DataGeneric as well as some attribute (_status for example)
- support 3w trafo (as modeled in pandapower)
- improve speed by not performing internal checks 
  (keep check for boundaries and all for python API instead) [see `TODO DEBUG MODE` in c++ code]
- improve speed
- code parrallelism directly in the `Computer` and `SecurityAnalysisCPP` classes
- a mode to do both `Computer` and `SecurityAnalysisCPP`
- use the "multi slack hack" (see issue #50) for SecurityAnalysis or Computer for example
- code `helm` powerflow method
- interface with gridpack (to enforce q limits for example)
- maybe have a look at suitesparse "sliplu" tools ?
- easier building (get rid of the "make" part)

TODO: speed directly update the pv, pq, Sbus and Ybus part when updating the elements
      (less error prone and faster to recompute). Then what is passed to the solver 
      is "only" a "mask" of these data when null
TODO: https://github.com/haranjackson/NewtonKrylov for another type of algorithm ?
TODO HVDC in Jacobian (see pandapower)
TODO: in ContingencyAnalysisCpp: add back the `if(!ac_solver_used)` inside the  `remove_from_Ybus`
      in order to perform the "invertibility" check
TODO: in `main.cpp` check the returned policy of pybind11 and also the `py::call_guard<py::gil_scoped_release>()` stuff
TODO: a cpp class that is able to compute (DC powerflow) ContingencyAnalysis and TimeSeries using PTDF and LODF
TODO: integration test with pandapower (see `pandapower/contingency/contingency.py` and import `lightsim2grid_installed` and check it's True)

[0.12.2] 2026-02-xx
----------------------
- [FIXED] an issue with shunt buses (was set to 1 even if they were disconnected)
- [FIXED] a warning when applying actions on generator votlage setpoints (due to NaN)
- [IMPROVED] add a test to make sure generator types are available if using
  `dist_slack_non_renew` information.

[0.12.1]  2026-01-09
---------------------
- [FIXED] phase shift transformers are now properly modeled
  for both pandapower (new in this version) and pypowsybl (already
  the case in previous version)
- [FIXED] a performance issue for all "XXXSingleSlack" (*eg* KLUSingleSlack) algorithm (filling of the initial
  Jacobian matrix was extremly slow due to the massive 'insert' of data in the eigen sparse matrix instead 
  of relying on the "setFromTriplets" method)
- [FIXED] an normal "exception" was not catched in the `close()` method of LightSimBackend in case the
  backend was closed before any grid was loaded.
- [ADDED] possibility to pickle independantly all part of the grid (*eg* gridmodel.get_lines()
  can be pickled independantly from anything else) **NB** pickling and un-pickling 
  lightsim2grid objects can only be used for the same lightsim2grid version.
- [ADDED] the `init_from_n_powerflow` property for ContingencyAnalysis (and 
  ContingencyAnalysisCPP). It allows to chose if the computation of the contingencies
  are initialized with the complex voltages resulting of the powerflow in N 
  (`init_from_n_powerflow=True`) or if they are initialized from the 
  given input vector (`init_from_n_powerflow=False`). Defaults to `False`

[0.12.0] 2026-01-06
--------------------
- [BREAKING] for better consistency, and following pypowsybl convention, trafo and lines "side"
  are now called "1" and "2" instead of "hv" / "lv" (for trafo) or "or" / "ex" for powerlines.
  For example, what used to be accessible with `gridmodel.change_bus_powerline_or(...)` is now called
  `gridmodel.change_bus1_powerline()`. This affects powerlines, transformers and dc powerlines but also
  "LineInfo", "TrafoInfo" and "DCPowerlineInfo". See below for a (should-be exhaustive) list of changes.
- [BREAKING] the `init_pp_backend` public attribute is now private (called now `_init_pp_backend`) and 
  optional, meaning it's `None` when the grid is initialized from pypowsybl (for example).
- [FIXED] some issues with the "load_grid_from_pypowsybl" function (and making sure the graph of the structure
  of the lightsim2grid gridmodel matches the one of the pypowsybl grid).
- [FIXED] an issue with the handling of the slack due to a not correct implementation
  of `update_slack_weights_by_id` cpp side (previous slacks were not removed, slacks get_slack_weights
  were prop to target_p which caused issues when all slacks had targetp==0.)
- [FIXED] an issue with serialization / de serialization caused by an error in serializing the solver types.
- [ADDED] in all "xxxInfo"  (*eg* "LoadInfo") information about subtation and position in the topology
  vector, with the `sub_id` / `pos_topo_vect` (for `LoadInfo`, `SGenInfo`, `GenInfo`, `StorageInfo`, `ShuntInfo`)
  and `sub1_id` / `sub2_id` / `pos1_topo_vect` / `pos2_topo_vect` (for `LineInfo`, `TrafoInfo` and `DCLineInfo`)
- [ADDED] possibility to load the pypowsybl grid with extra key-words arguments (by using `pypowsybl_load_kwargs` in the 
  `loader_kwargs` of LightSimBackend)
- [ADDED] possibility to initialize LightSimBackend with an already loaded grid (by using the `grid` key of the
  `loader_kwargs` of LightSimBackend when loading it with pypowsybl)
- [ADDED] possibility to change the ratio (`rho`) of transformers (`gridmodel.change_ratio_trafo(trafo_id, new_rho)`)
- [ADDED] possibility to change the phase shift (`alpha`) of transformers (`gridmodel.change_shift_trafo(trafo_id, new_alpha)`)
- [ADDED] more consistency checkings to avoid "negative buses" cpp side.
- [ADDED] information about the coefficients assigned on the Ybus matrix for `LineInfo` and `TrafoInfo`: `yac_11`, `yac_12`, 
  `yac_21`, `yac_22`, `ydc_11`, `ydc_12`, `ydc_21`, and `ydc_22`
- [ADDED] possibility to have a powerline / transformer connected on only one side (for DC and Newton-Raphson algorithm, not implemented
  for fast-decoupled yet). This means that powerlines / transformers have 3 statuses: one for each side and one "global".
- [ADDED] possibility to choose the way lightsim2grid will internally treat the powerlines status:

  - `ignore_status_global` (`gridmodel.set_ignore_status_global(True)` or `gridmodel.set_ignore_status_global(False)`). If set
     to `False` (default) the the "global" status is synch with all the others. For example you can deactivate both sides of a line by 
     deactivating `status_global` and conversely if you deactivate both side of a given line, then its "status_global" is set 
     to "activated". If `ignore_status_global` is set to `True` then `global_status` is not updated at all and ignored (*NB* in 
     this case, calling `gridmodel.deactivate_line(...)` will deactivate both sides but not the "global status")
  - `synch_status_both_side` (`self.model.set_synch_status_both_side(XXX)`). If set to `True` (default) then the status of each 
    side of any given line will be synched. Meaning that if an action disconnects one side, it will also disconnect the "status_global"
    and the other side (*NB* in this case if the same actions both connects one side and disconnect another, then the outcome is "undefined").
    If `synch_status_both_side` is `False` then each side of the powerline is independant from the other (which can lead to powerline / 
    transformer being connected at only one side).
  - The complete bahviour is tested in `tests/test_line_disco_one_side.py`. Feel free to have a look if you need more information.

- [IMPROVED] Eigen to version 5.0.1 (2025/11/11)
- [IMPROVED] rename all ".h" file to ".hpp" for cpp headers (cpp side). 
- [IMPROVED] consistency of "bus labelling" cpp side (implement compile time check to prevent accidental conversion from 
  `LocalBusId`, `GlobalBusId` / `GridModelBusId` and / or `SolverBusId`)

Table to upgrade the names:

=============================    =======================  =======================
Class Name                       Old Attribute Name       New Attribute Name
=============================    =======================  =======================
TrafoContainer                   get_bus_from             get_bus_id_side_1
TrafoContainer                   get_bus_to               get_bus_id_side_2
TrafoInfo                        bus_hv_id                bus1_id
TrafoInfo                        bus_lv_id                bus2_id
TrafoInfo                        connected                connected_global
TrafoInfo                        h_pu                     h1_pu or h2_pu
TrafoInfo                        is_tap_hv_side           is_tap_side_1
TrafoInfo                        res_p_hv_mw              res_p1_mw
TrafoInfo                        res_q_hv_mvar            res_q1_mvar      
TrafoInfo                        res_v_hv_kv              res_v1_kv    
TrafoInfo                        res_a_hv_ka              res_a1_ka         
TrafoInfo                        res_p_lv_mw              res_p2_mw      
TrafoInfo                        res_q_lv_mvar            res_q2_mvar    
TrafoInfo                        res_v_lv_kv              res_v2_kv     
TrafoInfo                        res_a_lv_ka              res_a2_ka     
TrafoInfo                        res_theta_hv_deg         res_theta1_deg    
TrafoInfo                        res_theta_lv_deg         res_theta2_deg
LineContainer                    get_bus_from             get_bus_id_side_1
LineContainer                    get_bus_to               get_bus_id_side_2
LineInfo                         connected                connected_global
LineInfo                         bus_or_id                bus1_id
LineInfo                         bus_ex_id                bus2_id
LineInfo                         h_pu                     removed
LineInfo                         h_or_pu                  h1_pu
LineInfo                         h_ex_pu                  h2_pu
LineInfo                         res_p_or_mw              res_p1_mw
LineInfo                         res_q_or_mvar            res_q1_mvar      
LineInfo                         res_v_or_kv              res_v1_kv    
LineInfo                         res_a_or_ka              res_a1_ka         
LineInfo                         res_p_ex_mw              res_p2_mw      
LineInfo                         res_q_ex_mvar            res_q2_mvar    
LineInfo                         res_v_ex_kv              res_v2_kv     
LineInfo                         res_a_ex_ka              res_a2_ka     
LineInfo                         res_theta_or_deg         res_theta1_deg    
LineInfo                         res_theta_ex_deg         res_theta2_deg
DCLineContainer                  get_bus_from             get_bus_id_side_1
DCLineContainer                  get_bus_to               get_bus_id_side_2
DCLineInfo                       connected                connected_global
DCLineInfo                       bus_or_id                bus1_id
DCLineInfo                       bus_ex_id                bus2_id
DCLineInfo                       target_vm_or_pu          target_vm1_pu         
DCLineInfo                       target_vm_ex_pu          target_vm2_pu 
DCLineInfo                       gen_or                   gen1         
DCLineInfo                       gen_ex                   gen2
DCLineInfo                       res_p_or_mw              res_p1_mw
DCLineInfo                       res_q_or_mvar            res_q1_mvar      
DCLineInfo                       res_v_or_kv              res_v1_kv          
DCLineInfo                       res_p_ex_mw              res_p2_mw      
DCLineInfo                       res_q_ex_mvar            res_q2_mvar    
DCLineInfo                       res_v_ex_kv              res_v2_kv      
DCLineInfo                       res_theta_or_deg         res_theta1_deg      
DCLineInfo                       res_theta_ex_deg         res_theta2_deg      
=============================    =======================  =======================


[0.11.0] 2025-12-09
----------------------
- [DEPRECATED] python 3.7 builds will no longer be available
- [FIXED] a bug in the import of the grid from pypowsybl
- [FIXED] a bug with phase shifters in case the tap was on the low voltage side
- [FIXED] a bug with active shunt values (wrong sign in the cpp part) and wrong Ybus diagonal coeff
- [FIXED] a bug in DC computation with shunt active values (wrong sign)
- [FIXED] a bug in DC computation with some phase shifters (when tap was not tagged on correct side)
- [FIXED] a bug in FDPF when phase tap changer was not on high voltage side
- [FIXED] a lots of bug in the conversion of pypowsybl grid when using 
  "old" pypowysbl versions.
- [ADDED] compatibility with python 3.14 and python 3.14 build
- [ADDED] compatibility with pandapower >= 3 version when loading a grid 
  (pandapower changed the way it initilizes the transformers model parameters)
- [ADDED] more kwargs arguments are possible in the LightSimBackend `loader_kwargs`
- [ADDED] name of the substations are now read from the grid when initializing from 
  pypowsybl.
- [ADDED] support for multiple slack when reading a grid from pypowsybl.
- [IMPROVED] the way to initialize the transformers from pypowsybl
- [IMPROVED] possibility to load grid with phase shifters from pypowsybl
- [IMPROVED] function to initialize the grid from pypowsybl has now a 
  basic documentation.

[0.10.3] 2025-04-28
----------------------
- [FIXED] remove deprecated use of numpy<2 function in LightSimBackend
- [FIXED] the "copy.deepcopy" of a lightsim2grid backend does not crash anymore 
  (see issues #36 and #97)
- [FIXED] an issue that could lead to a "segfault" if the "runpf" method of
  `LightSimBackend` was called first with is_dc=True and then with is_dc=False
- [IMPROVED] compat with grid2op 1.11.0
- [IMPROVED] now test proper compilation on clang 20 (was limited to clang 18 before)

[0.10.2] 2025-03-07
----------------------
- [FIXED] an error when changing of bus one of the slack (did not trigger the 
  recompute of pv bus ids)
- [FIXED] an issue when turning off a generator: it was still declared as "slack"
  if it was one.
- [FIXED] could not disconnect a generator when it was a slack bus
- [FIXED] voltage was -1 instead of 0 for disconnected elements (load, generator, storage units etc.)
- [ADDED] an option in `LightSimBackend` automatically disconnect load and generators
  if they are not in the main connected component.
- [IMPROVED] refactoring of the c++ side container element to reduce
  code (for "one end" elements such as loads, generators, static generators and shunts)

[0.10.1] 2025-01-04
----------------------------
- [FIXED] some timings on the benchmarks were not measured at the right time
- [ADDED] more benchmarks especially for DC powerflow
- [ADDED] a `dcpf` function that can replace the pandapower `dcpf` interal function
- [IMPROVED] benchmark on the documentation
  (clarity of what is done)
- [IMPROVED] consistency of the names and measured times accross the different benchmarks

[0.10.0] 2024-12-17
-------------------
- [BREAKING] disconnected storage now raises errors if some power is produced / absorbed, when using legacy grid2op version,
  you can retrieve the previous behaviour by initializing the `LightSimBackend` with
  `backend = LightSimBackend(..., stop_if_storage_disco=False)`
- [BREAKING] with the new `detachment_is_allowed` feature in grid2op, the kwargs `stop_if_load_disco`,
  `stop_if_gen_disco` (and `stop_if_storage_disco`) are now optional. They are set up from the 
  call to `grid2op.make(...)` and are erased by the `allow_detachment` kwargs. In other words,
  you don't need to set `stop_if_load_disco`, `stop_if_gen_disco` or `stop_if_storage_disco`. It is 
  automatically set by `grid2op.make(..., allow_detachment=XXX)` to have the correct bahaviour.
- [FIXED] an issue with the storage units (when asking it to produce / consume 
  but deactivating them with the same action the grid did not diverge)
- [IMPROVED] add the grid2op "detachement" support (loads and generators are allowed
  to be disconnected from the grid)
- [ADDED] a kwargs `stop_if_storage_disco` to control (in legacy grid2op version) the behaviour 
  of the backend when a storage unit is disconnected.

[0.9.2.post2] 2024-11-29
--------------------------
- [FIXED] The attribute `can_output_theta` (of base `Backend` class)
  was not set to `True` if using the pypowsybl loader.
- [FIXED] the github CI to work properly on many linux buit image
- [IMPROVED] build on python 3.13

[0.9.2.post1] 2024-11-28
--------------------------
- [FIXED] There is still a bug with the pypowsybl 1.8.1 version with the 
  tap changer ratio (unconsistency between what needs to be done and the 
  actual documentation). The fix is to set the const variable `PP_BUG_RATIO_TAP_CHANGER`
  to be at least 1.9.0, otherwise results are wrong.
  
[0.9.2] 2024-10-18
--------------------------
- [ADDED] support loading a grid when everything is NOT on the same bus
  (`topo_vect` used to be wrong in this case). This is especially usefull
  for grid loaded with `pypowsybl`
- [ADDED] a file benchmarking the timings for running powerflow on different
  grid sizes.
- [UPDATED] urls to match the new repo location
- [UPDATED] urls to match new grid2op location

[0.9.1] 2024-09-30
--------------------------
- [FIXED] a bug due to wrong type (in a numpy array) for the element name which lead in turn 
  to a fail assertion (equality between two numpy arrays returning a bool and not an array)
- [FIXED] a bug when init a grid from pypowsybl: the wrong value was used for trafos `h` (double)
- [FIXED] a bug when init a grid from pypowsybl: wrong values for `_ls_to_orig` and `_orig_to_ls`
  was set (and later used)
- [FIXED] yet another bug when init a grid from pypowsybl: the voltage in kV (not in pu)
  could be set due to "wrong" labelling of the bus ids
- [FIXED] yet another bug when init a grid from pypowsybl: the ratio of the transformers 
  sent in lightsim2grid did not take into account the "`rated_u1` `rated_u2`" on both side 
  (only used on one side)
- [FIXED] yet another bug when init a grid from pypowsybl: the ratio of the transformers 
  sent in lightsim2grid did not take into account the ratio in the  `pypow_net.get_ratio_tap_changers()`
- [ADDED] a method for the `ContingencyAnalysisCPP` class that returns, for all contingencies
  in the contingency list, which will be simulated and which causes the grid to be disconnected.
- [ADDED] it is now possible to use "one substation" (voltage level) pypowsybl side is
  one substation in lightsim2grid.
- [IMPROVED] removing a weird `1j * h_` when initializing powerlines and transformers. This was 
  part of a pandapower "hack" which is not present anymore (see 
  https://github.com/Grid2Op/lightsim2grid/issues/88#issue-2443299039)

[0.9.0] 2024-07-29
--------------------------
- [BREAKING] installing pandapower lightsim2grid does not require anymore to install
  pandapower (you can initialize `GridModel` with pypowsybl or pandapower if you want). To make it both
  cleaner and clearer the function `lightsim2grid.gridmodel.init` has been removed.
  Please use `lightsim2grid.gridmodel.init_from_pandapower` or 
  `lightsim2grid.gridmodel.init_from_pypowsybl` from now on
- [BREAKING] the previous `gridmodel.get_ptdf()` function was wrongly labelled with the
  "solver" bus id and not the `gridmodel` bus id which could cause issue when it was computed 
  on some grid configuration. It has now been fixed (so the `gridmodel.get_ptdf` returns the
  proper things). If you want the previous behaviour, you need to use `gridmodel.get_ptdf_solver()`
- [BREAKING] similarly, `gridmodel.get_Ybus()`, `gridmodel.get_dcYbus()`, `gridmodel.get_Sbus()`
  and `gridmodel.get_dcSbus()` now return things in the `gridmodel` bus ordering. For the previous
  behaviour you can use `gridmodel.get_Ybus_solver()`, `gridmodel.get_dcYbus_solver()`,
  `gridmodel.get_Sbus_solver()` and `gridmodel.get_dcSbus_solver()`
- [BREAKING] the more rational logic above also extends to all the functions listed in the 
  table below:

===============================    ===================================================
Function with behaviour change      Name of the new function having the same behaviour
===============================    ===================================================
gridmodel.get_ptdf()                gridmodel.get_ptdf_solver()
gridmodel.get_Ybus()                gridmodel.get_Ybus_solver()
gridmodel.get_dcYbus()              gridmodel.get_dcYbus_solver()
gridmodel.get_Sbus()                gridmodel.get_Sbus_solver()
gridmodel.get_dcSbus()              gridmodel.get_dcSbus_solver()
gridmodel.get_pv()                  gridmodel.get_pv_solver()
gridmodel.get_pq()                  gridmodel.get_pq_solver()
gridmodel.get_slack_ids()           gridmodel.get_slack_ids_solver()
gridmodel.get_slack_ids_dc()        gridmodel.get_slack_ids_dc_solver()
gridmodel.get_slack_weights()       gridmodel.get_slack_weights_solver()
gridmodel.get_V()                   gridmodel.get_V_solver()
gridmodel.get_Va()                  gridmodel.get_Va_solver()
gridmodel.get_Vm()                  gridmodel.get_Vm_solver()
gridmodel.get_J()                   gridmodel.get_J_solver()
gridmodel.get_Bf()                  gridmodel.get_Bf_solver()
===============================    ==================================================

- [FIXED] the `change_solver` in the `ContingencyAnalysis` did not work correctly.
  More specifically the solver type used might not be correct if changed which could 
  lead to wrong Ybus being passed to the solver.
- [FIXED] some compatibility mode with python `3.7`
- [FIXED] a bug when "turned off" generator were not PV (slack was 
  "turned off" when its target P was 0. But still the slack so ends up producing something...)
- [FIXED] (consistency with pandapower) when an intial powerflow is run
  to initialized an AC powerflow, the initial voltages are 1 pu (and 
  not `gridmodel.get_init_vm_pu()` as previously).
- [FIXED] `gridmodel.get_ptdf()` now have the 
  normal "gridmodel" bus id representation and not the "solver" bus ordering.
- [FIXED] `gridmodel.get_lodf()` issue wrong results in case of some
  topological modification
- [FIXED] calls to methods such as `gridmodel.get_pv` or `gridmodel.get_V` 
  or `gridmodel.get_Ybus` could lead to severe crashes (segmentation fault)
  on some (rare) cases. Now an exceptions should be thrown in these cases.
- [FIXED] basic backward compatibility is ensured and tested for legacy grid2op >= 0.9.1.post1
  Not all features are tested and only 1.x versions are tested 
  (ie 1.1 or 1.2 but not 1.2.1, 1.2.2, 1.2.3 etc.) and only for python 3.11
- [FIXED] a bug when using `LightSimBackend` with some old (but not too old) grid2op
  versions.
- [FIXED] various compatibility bugs when using old grid2op versions.
- [ADDED] it is now possible to deactivate the support for shunts by 
  subclassing the LightSimBackend class and setting the `shunts_data_available`
  to `False`
- [IMPROVED] in the `ContingencyAnalysis` class, the underlying cpp model will now
  perform an initial powerflow.
- [IMPROVED] distributed wheels are now compiled (whenever possible) with numpy 2. 
  This makes them compatible with both numpy 1.x.y and numpy 2.z.t versions.
- [IMPROVED] tests are now performed when lightsim2grid is compiled with 
  the latest clang (18) and gcc (14) versions on the CI using python 3.11

[0.8.2.post1] 2024-04-xx
--------------------------
- [FIXED] a "forward compatibility" issue with grid2op 1.10.2
  (due to wrong usage of some internal classes when loading a pandapower grid)

[0.8.2] 2024-04-22
--------------------
- [FIXED] CI was broken after migration to artifact v4, set it back to v3 
  (and make the names of the folder clearer)
- [FIXED] CI when using latest pandapower version (2.14) which broke some previous tests
- [ADDED] the computation of the LODF (line outage distribution factor) in 
  lightsim2grid
- [ADDED] some convenience functions to retrieve in a vectorized way the 
  buses to which each elements of a given container is connected 
  (*eg* `gridmodel.get_lines().get_bus_from()`)
- [ADDED] more binaries (windows `arm64` and macos `arm64`)
- [IMPROVED] remove some compilation warnings for clang
- [IMPROVED] possibility to specify generator used as slack by its name when initializing
  from `pypowsybl`.
- [IMPROVED] removing some warnings when grid2op is not installed
  (it should not raise any warning as lightsim2grid does not require grid2op)

[0.8.1] 2024-03-26
--------------------
- [FIXED] a bug with shunts when `nb_busbar_per_sub` >= 2
- [FIXED] some bugs preventing backward compatibility
- [FIXED] an issue in the computation of gen_q when intialized with pypowsybl
  (some overflow cpp side leading to infinite number in gen_q)
- [FIXED] a bug in the "containers" cpp side (wrong bus was assigned)
  when elements was disconnected, which lead to wrong computations for 
  time series or contingency analysis.
- [FIXED] another bug in ContingencyAnalysis (cpp side) leading to wrong computation
  when a powerline was disconnected
- [FIXED] some broken imports when grid2op was not installed
- [FIXED] missing "typing_extension" as required when installation
- [ADDED] some information of compilation directly in the cpp module
- [ADDED] some information of compilation available in the python `compilation_options`
  module python side
- [ADDED] some convenient methods for `ContingencyAnalysis` python side (most 
  notably the possibility to initialize it from a `LightSimBackend` and to
  change the topology of the grid)
- [ADDED] a "reward" module in lightsim2grid with custom reward
  based on lightsim2grid.
- [ADDED] a class `N1ContingencyReward` that can leverage lightsim2grid to 
  assess the number of safe / unsafe N-1.
- [IMPROVED] time measurments in python and c++
- [IMPROVED] now test lightsim2grid with oldest grid2op version
- [IMPROVED] speed, by accelerating the reading back of the data (now read only once and then
  pointers are re used)
- [IMPROVED] c++ side avoid allocating memory (which allow to gain speed python side too)
- [IMPROVED] type hinting in `LightSimBackend` for all 'public' methods (most 
  notably the one used by grid2op)
- [IMPROVED] now the benchmarks are more verbose (detailing some compilation options)

[0.8.0] 2024-03-18
--------------------
- [BREAKING] now able to retrieve `dcSbus` with a dedicated method (and not with the old `get_Sbus`).
  If you previously used `gridmodel.get_Sbus()` to retrieve the Sbus used for DC powerflow, please use
  `gridmodel.get_dcSbus()` instead.
- [DEPRECATED] in the cpp class: the old `SecurityAnalysisCPP` has been renamed `ContingencyAnalysisCPP`
  (you should not import it, but it you do you can `from lightsim2grid.securityAnalysis import ContingencyAnalysisCPP` now)
- [DEPRECATED] in the cpp class: the old `Computers` has been renamed `TimeSerieCPP`
  (you should not import it, but it you do you can `from lightsim2grid.time_serie import TimeSerieCPP` now)
- [FIXED] now voltage is properly set to 0. when shunts are disconnected
- [FIXED] now voltage is properly set to 0. when storage units are disconnected
- [FIXED] a bug where non connected grid were not spotted in DC
- [FIXED] a bug when trying to set the slack for a non existing genererator
- [FIXED] a bug in init from pypowsybl when some object were disconnected. It raises
  an error (because they are not connected to a bus): now this function properly handles
  these cases.
- [FIXED] a bug leading to not propagate correctly the "compute_results" flag when the 
  environment was copied (for example)
- [FIXED] a bug where copying a lightsim2grid `GridModel` did not fully copy it
- [FIXED] a bug in the "topo_vect" comprehension cpp side (sometimes some buses 
  might not be activated / deactivated correctly)
- [FIXED] a bug when reading a grid initialize from pypowsybl (trafo names where put in place 
  of shunt names)
- [FIXED] read the docs was broken
- [FIXED] a bug when reading a grid from pandapower for multiple slacks when slack 
  are given by the "ext_grid" information.
- [FIXED] a bug in "gridmodel.assign_slack_to_most_connected()" that could throw an error if a 
  generator with "target_p" == 0. was connected to the most connected bus on the grid
- [FIXED] backward compat with "future" grid2op version with a 
  better way to copy `LightSimBackend`
- [ADDED] sets of methods to extract the main component of a grid and perform powerflow only on this
  one.
- [ADDED] possibility to set / retrieve the names of each elements of the grid.
- [ADDED] embed in the generator models the "non pv" behaviour. (TODO need to be able to change Q from python side)
- [ADDED] computation of PTPF (Power Transfer Distribution Factor) is now possible
- [ADDED] (not tested) support for more than 2 busbars per substation
- [ADDED] a timer to get the time spent in the gridmodel for the powerflow (env.backend.timer_gridmodel_xx_pf)
  which also include the time 
- [ADDED] support for more than 2 busbars per substation (requires grid2op >= 1.10.0)
- [ADDED] possibility to retrieve the bus id of the original iidm when initializing from pypowsybl 
  (`return_sub_id` kwargs). This is a "beta" feature and will be adressed in a better way
  in a near future.
- [ADDED] possibility to continue the grid2op 'step' when the solver converges but a load or a 
  generator is disconnected from the grid.
- [IMPROVED] now performing the new grid2op `create_test_suite` 
- [IMPROVED] now lightsim2grid properly throw `BackendError`
- [IMPROVED] clean ce cpp side by refactoring: making clearer the difference (linear) solver
  vs powerflow algorithm and move same type of files in the same directory. This change
  does not really affect python side at the moment (but will in future versions)
- [IMPROVED] CI to test on gcc 13 and clang 18 (latest versions to date)
- [IMPROVED] computation speed: grid is not read another time in some cases.
  For example, if load and generators do not change, then Sbus is not
  recomputed. Likewise, if the topology does not change, then the Ybus 
  is not recomputed either see https://github.com/Grid2Op/lightsim2grid/issues/72

[0.7.5.post1] 2024-03-14
-------------------------
- [FIXED] backward compat with "future" grid2op version with a 
  better way to copy `LightSimBackend`
  
[0.7.5] 2023-10-05
--------------------
- [FIXED] a bug in DC powerflow when asking for computation time: it was not reset to 0. when
  multiple powerflows used the same solver
- [FIXED] a bug in AC and DC powerflow when shunts had active values
- [ADDED] possibility to initialize a powergrid based on pypowsybl 
  see https://github.com/Grid2Op/lightsim2grid/issues/53
- [ADDED] some more algorithm to perform powerflow: Fast Decoupled Powerflow (in BX and XB variant)
  see https://github.com/Grid2Op/lightsim2grid/issues/63
- [ADDED] build lightsim2grid for python 3.12
- [ADDED] support for non distributed slack but multiple slack buses
  see https://github.com/Grid2Op/lightsim2grid/issues/50 (ONLY FOR AC powerflow)
- [IMPROVED] now shipping `src` and `eigen` directory in the source of 
  lightsim2grid to allow their installation if wheels are not provided.
- [IMPROVED] in the underlying cpp GridModel powerlines can now have 2
  different values for the `h` parameters (`h_or` and `h_ex`).
- [IMPROVED] now lightsim2grid is able to load a pandapower network with non
  contiguous non starting at 0 bus index

[0.7.3/4] 2023-08-24
--------------------
- [FIXED] a bug where, when you disconnect a load (or gen), the next action cannot be performed
  if it modifies the load (or gen), because you "cannot change the value of a disconnected load (or gen)"
- [FIXED] read-the-docs template is not compatible with latest sphinx version (7.0.0)
  see https://github.com/readthedocs/sphinx_rtd_theme/issues/1463
- [IMPROVED] initialize the underlying "PandaPowerBackend" without numba
- [IMPROVED] grid2op import to be more compliant with renaming of uppercased file names
- [IMPROVED] decoupling of the PandapowerBackend class and the class "internally" used by LightSimBackend
  when loading the grid. This caused some issue, *eg* https://github.com/Grid2Op/grid2op/issues/508

[0.7.2] 2023-06-06
--------------------
- [FIXED] a bug in the `init` function that caused issue when importing a grid with multiple slack
  on some cases
- [FIXED] some bugs in the "SecurityAnalysis" and "TimeSerie" modules especially in DC mode.
- [FIXED] a bug in the DC comptuation: some "divergence" were not catched
- [FIXED] a bug in the "Computer" (cpp) class where the intial voltage could lead to generator not
  participating correctly to the voltage regulation (wrong output voltage level).
- [FIXED] a bug in the "set_bus" of shunt (wrong bus was assigned cpp side)
- [FIXED] an issue when slack bus is added from ext grid (wrong active power value - sign issue)
- [ADDED] support for the CKTSO linear solver (on linux), which is slightly faster than SparseLU, KLU and NICSLU
  (this requires a compilation from source)
- [ADDED] support for distributed slack bus in `LightSimBackend`
- [ADDED] support for "generator with p=0. do not participate in voltage regulation" in `LightSimBackend`
- [ADDED] support for the DC computation for "SecurityAnalysis" and "TimeSerie" modules
- [ADDED] support for DC powerline (in lightsim, they are still not handled in grid2op)
- [IMPROVED] now that multiple slacks is fully supported, the warnings when importing a grid with multiple slacks
  are irrelevant. They have been removed.
- [IMPROVED] the documentation on the "sovlers" part
- [IMPROVED] move the "how to compile" section of the readme in the documentation
- [IMPROVED] `SuiteSparse` is upgraded to version 5.13 (issue with build system based on cmake and BLAS for SuiteSparse >= 6.0)
- [IMPROVED] upgrade to eigen `3.4.0` (stable release)

[0.7.1] 2023-01-11
---------------------
- [BREAKING] drop support for numpy version < 1.20 (to be consistent with grid2op)
- [FIXED] a compatibility issue with grid2op 1.7.2 (missing another backend attribute
  when the environment is copied) see https://github.com/rte-france/Grid2Op/issues/360
- [FIXED] now an error if thrown if the bus indexes in the pandapower grid are not contiguous
  or do not start at 0 (thanks Roman Bolgaryn for spotting this issue)
- [ADDED] automatic build for python 3.11
- [ADDED] support for numpy >= 1.24 (some deprecation *eg** np.str and np.bool are removed)

[0.7.0.post1] 2022-06-20
-------------------------
- [FIXED] a compatibility issue with grid2op 1.7.1 (missing a backend attribute
  when environment is copied)

[0.7.0] 2022-05-30
---------------------
- [ADDED] improved time measurments
- [ADDED] Possibility to set, at creation time, the type of solver used, number
  of iterations and precisions with 
  `LightSimBackend(max_iter=..., tol=..., solver_type=...)`
- [IMPROVED] scripts to load the pandapower grid (json format)
- [IMPROVED] update the automatic tests on more recent compilers.

[0.6.1.post2] 2022-02-08
-------------------------
- [FIXED] add support for python 3.10 now that scipy does (and add proper tests in CI)

[0.6.1.post1] 2022-02-02
-------------------------
- [FIXED] support for python3.7 (and add proper tests in CI)

[0.6.1] 2022-02-01
--------------------
- [BREAKING] the behaviour of the `newton_pf` function is not 
  consistent with pandapower default concerning distributed slack.
- [FIXED] an issue in the distributed slack case spotted by pandapower team 
  thanks to them (see https://github.com/e2nIEE/pandapower/pull/1455)
- [IMPROVED] lightsim2grid will now use the single slack algorithm if the 
  grids counts only one slack bus (performance increase)

[0.6.0] 2021-12-17
-------------------
- [BREAKING] change the interface of the `newton_pf` function to reflect pandapower change in their
  latest version (arguments `ref` has been added). You can still use the old `newton_pf` function, with the
  old signature by importing `newtonpf_old` instead or explicitly importing the new one by importing `newtonpf_new`
- [BREAKING] `SecurityAnalysis` now also returns the active flows when calling `security_analysis.get_flows()`
- [BREAKING] change the file names (python side) to be compliant with pep 8. You can no longer
  do things like `from lightsim2grid.LightSimBackend import LightSimBackend` change it to
  `from lightsim2grid import LightSimBackend` (preferred method)
- [BREAKING] change the file names (python side) to be compliant with pep 8. You can no longer
  do things like `from lightsim2grid.initGridModel import init` change it to
  `from lightsim2grid.gridmodel import init` (preferred method) (same for `GridModel` class)
- [FIXED] a bug that lead to the wrong computation of the dc powerflow in case of `sn_mva != 1.` and phase shifters.
- [FIXED] bug preventing to use the NICSLU linear solver in the `GridModel`
- [FIXED] compilation warnings on clang (missing virtual destructor, unused variables, etc.)
- [FIXED] a bug in the `SecurityAnalysisCPP`: when it diverges for some contingencies, the others were not simulated properly.
- [FIXED] `LightSimBackend` now contains members for `shunts` and `***_theta` as it does for the other quantities. This improves the consistency, but most importantly
  fixes some bugs when used in earlier grid2op versions
- [ADDED] possibility to compute the active flows using the `BaseMultiplePower` 
- [ADDED] possibility to change linear solver used when performing a DC solver
- [ADDED] possibility to make powerflow with distributed slack bus (only for newton raphson at the moment)
- [ADDED] access (read only) to the element of a lightsim2grid grid with the `get_XXX` (*eg* `get_loads()`) methods (see documentation)
- [ADDED] direct access to the solver used in the grid model python side
- [ADDED] unittest in circleci.
- [ADDED] all kind of solvers based on different linear solvers (Eigen sparse LU, KLU or NICSLU) for Newton Raphson and
  DC approximation (9 solvers in total)
- [IMPROVED] use of `steady_clock` to retrieve the ellapse time c++ side
- [IMPROVED] refactoring of the c++ part to use template mecanism instead of inheritance for the
  Newton Raphson and DC solvers.
- [IMPROVED] `GridModel` now contains two different solvers, one for AC powerflow and one for DC powerflow.
- [IMPROVED] error message in the solver are now embedded in an Enum instead of being integers, for better readibility.
- [IMPROVED] error message when the powerflow diverge (error are read from c++ now)

[0.5.5] 2021-11-10
-------------------
- [ADDED] possibility to perform dc powerflow
- [ADDED] a class to compute flows on whole time series when the Ybus does not change (see `TimeSerie`)
- [ADDED] a class to compute flows on multiple contingencies, when Sbus does not change (see `SecurityAnalysis`).
- [IMPROVED] running speed of Newton Raphson solvers with better filling of sparse matrices
- [IMPROVED] upgrade to SuiteSparse `v5.10.1`
- [IMPROVED] upgrade to eigen `3.4.0` (stable release)
- [IMPROVED] clean the compilation warnings on microsoft windows (force the conversion from
  `Eigen::EigenBase<Derived>::Index` to `int` using `static_cast`)
- [IMPROVED] add the proper optimization flag for windows (`/O2` instead of `-03` on linux / macos)
- [IMPROVED] high performance gain when topology is not changed between steps (gain obtained by 
  reusing the previous Ybus)

[0.5.4] 2021-08-20
------------------
- [FIXED] a bug for static generator (wrong signed convention were used in some part of the c++ code). This has
  no impact at all for provided grid2op environments.
- [FIXED] An issue where the backend could get "stuck" in a wrong state because of the way the Vinit was computed (see
  `Issue 30 <https://github.com/Grid2Op/lightsim2grid/issues/30>`_)
- [ADDED] experimental support for the `NICSLU` linear solver (requires a proper license and library, see
  https://github.com/chenxm1986/nicslu for more information. Support does not include multi threaded at the moment).
- [IMPROVED] minor performance improvements for the solvers based on Newton Raphson (faster filling of the Jacobian
  matrix after the first iteration)

[0.5.3] 2021-08-11
-------------------
- [FIXED] minor issues in the benchmark (some time measurments were wrong)
- [ADDED] lightsim2grid package now can be distributed on pypi
- [ADDED] compilation of SuiteSparse using cmake
- [ADDED] compatibility with the KLU linear solver on windows based systems.
- [IMPROVED] the package should now be available on pypi

[0.5.2] 2021-07-26
-------------------
- [FIXED] `GridModel` now properly throw "out_of_range" exception when trying to change the bus of non existing
  elements
- [FIXED] wrong units were displayed for the iterators for lines and transformers.
- [ADDED] now able to retrieve the powerlines parameters python side.
- [IMPROVED] more explicit error messages when the building of the `Ybus` matrix fails.
- [IMPROVED] now the solver is not reset when using the `backend._grid.check_solution`
- [IMPROVED] upgrade SuiteSparse to version `v5.10.1`
- [IMPROVED] upgrade eigen to version `3.4-rc1`

[0.5.1] 2021-04-09
-------------------
- [FIXED] yet another compilation issue with clang (see
  `Issue 22 <https://github.com/Grid2Op/lightsim2grid/issues/22>`_)
- [ADDED] circleci to check compilation for gcc
- [ADDED] circleci to check compilation for clang
- [ADDED] circleci to check compilation for msvc
- [ADDED] function to read the voltage angle from the backend
- [ADDED] compatibility with grid2op 1.5.0 (up to an issue with the storage units)

[0.5.0] 2021-03-01
-------------------
- [FIXED] a compilation issue on macos
- [FIXED] a compilation issue on windows (missing import of vector in `DataConverter.h`)
- [FIXED] an import issue (with `lightsim2grid.SolverType`)
- [FIXED] a bug that lead to the wrong computation of the ratio of the trafo when the tap on hv side.
- [FIXED] wrong timing was measured in the "solver powerflow time" of pandapower in the benchmarks
- [FIXED] a broken handling of shunt modification (wrong bus was assigned)
- [FIXED] an issue in `LightSimBackend.copy` that prevent the copied environment from being reset.
- [FIXED] errors are now raised when pandapower grid cannot be converted in lightsim2grid (*eg.* when
  unsupported elements are present)
- [ADDED] a variant of the Gauss Seidel method which does the update in a "synchronous" fashion
- [ADDED] a function that, given a complex vector is able to check kicchoff's law violation.
- [ADDED] Support for phase shifter (modeled as trafo with an extra parameter `shift`)
- [ADDED] Experimental support for `sn_mva` pandapower parameter.
- [UPDATED] github issue template
- [IMPROVED] warnings are issued when some of the pandapowergrid attributes have been automatically replaced
  when converting to / from pandapower

[0.4.0] - 2020-10-26
---------------------
- [ADDED] the Gauss Seidel method for AC powerflow is now available
- [ADDED] possibility to change easily the solver types from python side

[0.3.0] - 2020-10-06
-------------------------
- [ADDED] Support for pickle for the lightsim Backend.
- [ADDED] LightSim should now be compatible with windows (implementation of a powerflow mode without
  using the SuiteSparse KLU linear solver but rather the Eigen SparseLU one)
- [ADDED] start of the documentation.

[0.2.4] - 2020-08-20
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
