Lightsim2grid 0.8.0 and grid2op 1.10.0
======================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-03-25 17:53  CET
- system: Linux 5.15.0-56-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-6820HQ CPU @ 2.70GHz
- python version: 3.10.13.final.0 (64 bit)
- numpy version: 1.23.5
- pandas version: 2.2.1
- pandapower version: 2.13.1
- grid2op version: 1.10.1
- lightsim2grid version: 0.8.1
- lightsim2grid extra information: 

	- klu_solver_available: True 
	- nicslu_solver_available: True 
	- cktso_solver_available: True 
	- compiled_march_native: True 
	- compiled_o3_optim: True 

====================  ======================  ===================================  ============================
case14_sandbox          grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
====================  ======================  ===================================  ============================
PP                                      46.3                               18.4                          6.57
GS                                     757                                  0.474                        0.378
GS synch                               769                                  0.445                        0.348
NR single (SLU)                        960                                  0.184                        0.0831
NR (SLU)                               952                                  0.189                        0.0819
NR single (KLU)                       1030                                  0.12                         0.0221
NR (KLU)                              1030                                  0.118                        0.0202
NR single (NICSLU *)                  1020                                  0.121                        0.022
NR (NICSLU *)                         1020                                  0.119                        0.02
NR single (CKTSO *)                   1020                                  0.119                        0.0211
NR (CKTSO *)                           989                                  0.121                        0.0192
FDPF XB (SLU)                         1010                                  0.13                         0.032
FDPF BX (SLU)                         1010                                  0.143                        0.0451
FDPF XB (KLU)                         1020                                  0.124                        0.0263
FDPF BX (KLU)                         1010                                  0.134                        0.0377
FDPF XB (NICSLU *)                    1010                                  0.126                        0.0267
FDPF BX (NICSLU *)                    1020                                  0.134                        0.0383
FDPF XB (CKTSO *)                     1010                                  0.125                        0.0268
FDPF BX (CKTSO *)                     1000                                  0.136                        0.0381
====================  ======================  ===================================  ============================

Differences:

============================  ==============  ==============  ================
case14_sandbox (1000 iter)      Δ aor (amps)    Δ gen_p (MW)    Δ gen_q (MVAr)
============================  ==============  ==============  ================
PP (ref)                            0               0                 0
GS                                  0.000122        7.63e-06          7.63e-06
GS synch                            0.000122        7.63e-06          7.63e-06
NR single (SLU)                     0.000122        7.63e-06          7.63e-06
NR (SLU)                            0.000122        7.63e-06          7.63e-06
NR single (KLU)                     0.000122        7.63e-06          7.63e-06
NR (KLU)                            0.000122        7.63e-06          7.63e-06
NR single (NICSLU *)                0.000122        7.63e-06          7.63e-06
NR (NICSLU *)                       0.000122        7.63e-06          7.63e-06
NR single (CKTSO *)                 0.000122        7.63e-06          7.63e-06
NR (CKTSO *)                        0.000122        7.63e-06          7.63e-06
FDPF XB (SLU)                       0.000122        7.63e-06          7.63e-06
FDPF BX (SLU)                       0.000122        7.63e-06          7.63e-06
FDPF XB (KLU)                       0.000122        7.63e-06          7.63e-06
FDPF BX (KLU)                       0.000122        7.63e-06          7.63e-06
FDPF XB (NICSLU *)                  0.000122        7.63e-06          7.63e-06
FDPF BX (NICSLU *)                  0.000122        7.63e-06          7.63e-06
FDPF XB (CKTSO *)                   0.000122        7.63e-06          7.63e-06
FDPF BX (CKTSO *)                   0.000122        7.63e-06          7.63e-06
============================  ==============  ==============  ================


l2rpn_neurips_2020_track2_small
---------------------------------

Configuration: 

- date: 2024-03-25 17:59  CET
- system: Linux 5.15.0-56-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-6820HQ CPU @ 2.70GHz
- python version: 3.10.13.final.0 (64 bit)
- numpy version: 1.23.5
- pandas version: 2.2.1
- pandapower version: 2.13.1
- grid2op version: 1.10.1
- lightsim2grid version: 0.8.1
- lightsim2grid extra information: 

	- klu_solver_available: True 
	- nicslu_solver_available: True 
	- cktso_solver_available: True 
	- compiled_march_native: True 
	- compiled_o3_optim: True 

=====================  ======================  ===================================  ============================
neurips_2020_track2      grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
=====================  ======================  ===================================  ============================
PP                                      41.5                                20.7                           8.6
GS                                       3.74                              266                           266
GS synch                                35.8                                26.9                          26.8
NR single (SLU)                        536                                   0.897                         0.767
NR (SLU)                               505                                   0.959                         0.818
NR single (KLU)                        811                                   0.268                         0.144
NR (KLU)                               820                                   0.256                         0.131
NR single (NICSLU *)                   813                                   0.259                         0.134
NR (NICSLU *)                          827                                   0.243                         0.118
NR single (CKTSO *)                    814                                   0.257                         0.131
NR (CKTSO *)                           829                                   0.24                          0.116
FDPF XB (SLU)                          762                                   0.352                         0.232
FDPF BX (SLU)                          749                                   0.373                         0.252
FDPF XB (KLU)                          786                                   0.307                         0.188
FDPF BX (KLU)                          776                                   0.327                         0.206
FDPF XB (NICSLU *)                     786                                   0.308                         0.188
FDPF BX (NICSLU *)                     771                                   0.324                         0.204
FDPF XB (CKTSO *)                      784                                   0.309                         0.19
FDPF BX (CKTSO *)                      773                                   0.329                         0.209
=====================  ======================  ===================================  ============================

Differences:

=================================  ==============  ==============  ================
neurips_2020_track2 (1000 iter)      Δ aor (amps)    Δ gen_p (MW)    Δ gen_q (MVAr)
=================================  ==============  ==============  ================
PP (ref)                                  0              0                 0
GS                                        6.1e-05        3.81e-06          1.53e-05
GS synch                                  6.1e-05        3.81e-06          1.53e-05
NR single (SLU)                           6.1e-05        0                 9.54e-07
NR (SLU)                                  6.1e-05        0                 9.54e-07
NR single (KLU)                           6.1e-05        0                 9.54e-07
NR (KLU)                                  6.1e-05        0                 9.54e-07
NR single (NICSLU *)                      6.1e-05        0                 9.54e-07
NR (NICSLU *)                             6.1e-05        0                 9.54e-07
NR single (CKTSO *)                       6.1e-05        0                 9.54e-07
NR (CKTSO *)                              6.1e-05        0                 9.54e-07
FDPF XB (SLU)                             6.1e-05        1.91e-06          1.53e-05
FDPF BX (SLU)                             6.1e-05        1.91e-06          7.63e-06
FDPF XB (KLU)                             6.1e-05        1.91e-06          1.53e-05
FDPF BX (KLU)                             6.1e-05        1.91e-06          7.63e-06
FDPF XB (NICSLU *)                        6.1e-05        1.91e-06          1.53e-05
FDPF BX (NICSLU *)                        6.1e-05        1.91e-06          7.63e-06
FDPF XB (CKTSO *)                         6.1e-05        1.91e-06          1.53e-05
FDPF BX (CKTSO *)                         6.1e-05        1.91e-06          7.63e-06
=================================  ==============  ==============  ================
