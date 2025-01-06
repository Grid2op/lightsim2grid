.. _benchmark-deep-dive:

Deep dive into the benchmarking of lightsim2grid
==================================================================

At various occasion, we benchmark lightsim2grid with other available solvers.

In most cases, we benchmark them using the `grid2op` package. In this section we briefly explain what happens and how to interpret the different figures and tables.

Grid2op in a nutshell
----------------------

This is a super basic overview of grid2op. Only the necessary to understand what is benchmarked is explained here. Please consult grid2op documentation at https://grid2op.readthedocs.io/ for more information.

When grid2op is used, a "grid2op environment" is run for a given number of steps (288 or 1000) usually with an "agent" that does nothing. From grid2op point of view this means:

1. loading the time series and initializing the solver to the initial state
2. for each of the steps (following steps repeated 288 or 1000 times):
   
   a. process the action of the agent (if the action is "do nothing" then nothing is really done here)
   b. retrieve the next value for all generators and load (from the "time series")
   c. "compile" everything into a "state" of the solver (state before solving for the Kirchhoff Current Laws - **KCL**)
   d. pass this information to the backend 
   e. inform the backend to run a powerflow (solve the Kirchhoff Current Laws - **KCL**)
   f. check for possible terminal condition (isolated loads or generators etc.)
   g. pass the results back to python and format them into a grid2op "observation"
3. stop the computation either if some terminal condition are reached (*eg* divergence of the powerflow) or if the time series have been all computed

.. note::
   This is a simplification of the possibility of grid2op. For example, by default grid2op can do multiple loops between steps 2e. and 2f. above.

   And once arrived at step 3. the process starts again at step 1. with a call to `env.reset`

Unless told otherwise, the `grid2op times`, `grid2op speed`, `step duration` and alike reports the average on the 288 (or 1000) points 2. above. 
This includes lots of computations (happening through grid2op) to read the data, compile them into a something that can be digested by the grid2op "backend"
and read back the results.

The time to compute the first and last steps (1. and 3. respectively) are never considered for these benchmarks.
Feel free to let us know (in the grid2op package) if this is of any interest to you.

All the other time reported are other divisions of the step 2e. above, which is detailed in the next section.


Grid2op Backend in a nutshell
-------------------------------

This section will detail the steps 2e. of the overview of the previous section, because it can, in turn, be decomposed into different steps.

The goal of this 2e. steps is to find the solution to the KCL for each buses of the grid. This is done thanks to a "solver". 
This solver is "wrapped" to grid2op with what is called a "grid2op backend". This is why there are 2 columns dedicated to this
"step 2e." in the page :ref:`benchmark-solver` for example, the column `grid2op 'backend.runpf' time (ms)` counts all the computation
happening in steps 2e. (and depending on the implementation - this is the case with lightsim2grid- even step 2f. and 2g.) whereas the 
column `solver powerflow time (ms)` only counts thet time spent in the "solver" for the AC powerflow.

What is happening in step 2e. (basically in the `backend.runpf` function) is:

1. before some basic initial steps (*eg* the connexity of the grid)
2. find some initial solution for the complex voltage at each bus 
   (this is done with the Direct Current approximation in pandapower and lightsim2grid)
3. start the resolution of the KCL in AC (for pandapower and lightsim2grid in most benchmarks)
4. derives all the others quantities from that (active and reactive flows on each line, 
   current flow on each line, reactive power absorbed / produced for each generators, voltage angle and magnitude 
   at each side of each powerline etc.)
5. pass all the data from the underlying data structure (*eg* "convert" from Eigen -a c++ library- vectors to numpy array if you use lightsim2grid)
6. check for possible terminal conditions that would stop the grid2op episode

Not all backends are forced to behave this way. For example, some backend might be initialized without first running a DC powerflow. 
Similarly, steps 5. and 6. above are not mandatory and can be done elsewhere in the code.

For lightsim2grid and pandapower (default implementation) all these steps are performed inside the `backend.runpf` function which makes them 
comparable to this regard.

At this stage, we know that `grid2op 'backend.runpf' time (ms)` corresponds to all the time spent in 2e including all 
the time to perform 2e1, 2e2, 2e3, 2e4, 2e5 and 2e6.

For pandapower and lightsim2grid, `solver powerflow time (ms)` does not report time spent on 2e3 (point 3. above: as you might 
remember this is a zoom into the 2e. steps of grid2op) but only a part of it (explained in the next, and last, section)

.. note:: 
   As of grid2op 1.11.0 (still in development when this page is being written), some check will disappear from pandapower and from lightsim2grid, especially 
   things checked at 2e6 (point 6 above)


Some vocabulary is still needed to explain some concepts of these benchmarks, especially the different solver used or the "recycling" term. For this, 
we need to dive deeper into the step 3 above (so diving deeper into 2e3 if we label things from the grid2op perspective).

Physical solver in a nutshell
------------------------------

All of the solvers (to solve these type of problems) that we are aware of actually have themselves two (or more) "layers".

There is the "external layer" that can be accessed and modified more or less easily. When these "physical solvers" perform powerflows (so the step 2e3, 
*ie* solving the KCL in AC), they will often : 

a. perform some check on the external layer / data model
b. convert this "external model" / "data model" into something that can be processed by an algorithm
c. run this algorithm
d. convert back the results into the "external layer" / "data model" so that user can easily access it


Now we can properly explain what is reported on the column `solver powerflow time (ms)` :

- for pandapower and lightsim2grid, it gives only the time spent in 2e3c (thus removing all the time spent in a, b and d). More precisely:

    - in lightsim2grid backend this is retrieved with `env.backend._grid.get_computation_time()` 
    - in pandapower backend it is obtained with `env.backend._grid["_ppc"]["et"]`
- for pypowsybl unfortunately, we do not have that much detail at hand. So the time reportd in `solver powerflow time (ms)` will 
  include all steps 2e3a, 2e3b, 2e3c and 2e3d.

In the pages :ref:`benchmark-solver` and :ref:`benchmark-grid-size` the concept of `recycling` is used without a proper definition. With this view
of the physical solver, we can start to explain it. A first type of "recycling" is to reuse previous data when the conversion between a. and b. happens. 
This can saves a lot of time.

For example, this can mean that:

- if the topology is not changed, then the `Ybus` matrix will not be recomputed
- if the injection is not changed, then the `Sbus` vector will not be recomputed
- memory will not be deallocated / reallocated between b -> c or between c -> d if possible (to save sometimes expensive system calls)
- some checks will not be done if the underlying data are not modified (skipping partially or totally step a)
- the algorithm (see next section) itself can "know" what part of its data can be reuse avoiding even further unnecessary computations
  (that is the "recycle property" can be also forwarded to the solver)
- etc. 

.. note::
   The split of a "physical solver" into these 4 steps holds for different "solvers", for example:
   
     - `pandapower`, where the "external layer" / "data model" consists of the dataframes reprensenting the grid
       (*eg* `pp_net.line`, `pp_net.load`, etc. when for the tables that can be modified or 
       `pp_net.res_line`, `pp_net.res_load`, etc. for read-only attribute)
     - `pypowsybl` where the "external layer" also consists of dataframes, accessible with `pp_grid.get_lines()` or `pp_grid.get_loads()`
       and the data within this model can be modified with `pp_grid.update_lines` or `pp_grid.update_loads`
     - `lightsim2grid` where the "gridmodel" data can be inspected with *eg* `gridmodel.get_lines()` or `gridmodel.get_loads()` And
       modified with `gridmodel.update_loads_p`, `gridmodel.update_loads_q` or `gridmodel.update_topology`


.. note::
   The stopping criteria of the algorithm might slightly differ depending on the "physical solver" used. We made sure that the 
   same stopping criteria is used for pandapower and lightsim2grid but it might differ for other solver.

Algorithm in a nutshell
---------------------------

This is the last step useful to understand the benchmarks performed (at time of writing) in lightsim2grid and it aims at 
understanding the last part of the `recycling` mecanism as well as questions like "why is there so much different rows in the solver benchmarks ?"

There are different algorithms to solve the AC KCL (the operation performed at step `2e3c`) among which:

- Gauss Seidel 
- Newton Raphson
- Fast Decoupled

When we benchmark lightsim2grid in the page :ref:`benchmark-solver` all 3 algorithms are tested, and for the page 
:ref:`benchmark-grid-sizebenchmark-grid-size` only Newton-Raphson algorithm is used (if you are interested in more, please
let us know, no problem at all).

Pandapower
++++++++++++

When pandapower is benchmarked, only the Newton-Raphson algorithm is used, we will not detail the exact implementatoin of pandapower. Its implementation
the python scipy package to perform the linear algebra operations needed.

Pypowsybl
++++++++++++

When pypowsybl is benchmarked, only the Newton-Raphson algorithm is used. It internally uses some java implementation relying on powsybl framework and
open-loadflow (for the default parameters) powerflow.

Let us know if you are interested with more detail and more algorithm (powsybl can do much more than what is exposed here).

Lightsim2grid
++++++++++++++

In the benchmarks, lightsim2grid counts the most reported algorithms. In this section we detail a "concisely" the bahviour all some of them. 

Gauss Seidel
~~~~~~~~~~~~~

Lightsim2grid comes with two different Gauss-Seidel algorithms. They are not very efficient for the kind of problem at hand, so 
we will not spend lot of time discussing them here.


Fast decoupled and Newton Raphson
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These two types of algorithms comes each with different variants (we will not enter into too their detail here):

- Newton-Raphson with a single slack bus
- Newton-Raphson with a distributed slack bus
- Fast Decoupled 'BX'
- Fast Decoupled 'XB'

And for each of these algorithm, at some point linear systems in the form of "Ax = b" (A and b known, solve for x, A being a sparse matrix) 
needs to be solved repeatedly. They can be decomposed in different steps (as always in this page
without entering into detail and with lots of simplifications):

1. perform some initial checks (to make sure data are consistent)
2. initialize the linear solver (require to allocate some memory, create some vectors, etc.)
3. repeat :
   1. check if maximum number of allowed iterations is reached (divergence) if so, stop
   2. check if the stopping criteria are met (convergence) if so, stop
   3. update a linear system based on the value of complex voltages at each buses
   4. solve the new linear system
   5. update the complex voltages at each buses

.. note::
   For Fast Decoupled method, actually two different types of update are made, this would be equivalent to having
   first do "3 and 4" and then do "3' and 4'" before finally updating the complex vector at step 5

.. note::
   The checks 1. and 2. might actually happen at the end of the loop (so after 5. in this case) depending on the algorithm
   but this does not change the message.

TODO recycling !!!  