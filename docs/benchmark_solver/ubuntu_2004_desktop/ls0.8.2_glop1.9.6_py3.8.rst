Lightsim2grid 0.8.2 and grid2op 1.9.6 (python 3.8)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 10:43  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.6
- lightsim2grid version: 0.8.2
- lightsim2grid extra information:
	- klu_solver_available: True
	- nicslu_solver_available: True
	- cktso_solver_available: True
	- compiled_march_native: True
	- compiled_o3_optim: True
	
====================  ======================  ===================================  ============================
case14_sandbox          grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
====================  ======================  ===================================  ============================
PP                                      32.9                               24                           16.4
GS                                     862                                  0.422                        0.331
GS synch                               854                                  0.431                        0.341
NR single (SLU)                       1100                                  0.16                         0.0661
NR (SLU)                              1110                                  0.161                        0.0668
NR single (KLU)                       1180                                  0.11                         0.019
NR (KLU)                              1170                                  0.111                        0.0177
NR single (NICSLU *)                  1150                                  0.112                        0.0193
NR (NICSLU *)                         1190                                  0.108                        0.0174
NR single (CKTSO *)                   1100                                  0.116                        0.0193
NR (CKTSO *)                          1090                                  0.117                        0.0181
FDPF XB (SLU)                         1170                                  0.118                        0.0283
FDPF BX (SLU)                         1150                                  0.13                         0.0405
FDPF XB (KLU)                         1170                                  0.114                        0.0235
FDPF BX (KLU)                         1160                                  0.123                        0.0333
FDPF XB (NICSLU *)                    1140                                  0.118                        0.025
FDPF BX (NICSLU *)                    1070                                  0.134                        0.0376
FDPF XB (CKTSO *)                     1080                                  0.125                        0.0263
FDPF BX (CKTSO *)                     1070                                  0.134                        0.037
====================  ======================  ===================================  ============================

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

- date: 2024-04-23 10:49  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.6
- lightsim2grid version: 0.8.2
- lightsim2grid extra information:
	- klu_solver_available: True
	- nicslu_solver_available: True
	- cktso_solver_available: True
	- compiled_march_native: True
	- compiled_o3_optim: True
	
=====================  ======================  ===================================  ============================
neurips_2020_track2      grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
=====================  ======================  ===================================  ============================
PP                                      15.7                                27.9                         19.1
GS                                       4.98                              200                          200
GS synch                                34.9                                27.7                         27.6
NR single (SLU)                        618                                   0.754                        0.632
NR (SLU)                               613                                   0.764                        0.642
NR single (KLU)                        911                                   0.238                        0.121
NR (KLU)                               918                                   0.228                        0.11
NR single (NICSLU *)                   966                                   0.22                         0.109
NR (NICSLU *)                          977                                   0.209                        0.0982
NR single (CKTSO *)                    968                                   0.216                        0.105
NR (CKTSO *)                           980                                   0.204                        0.0926
FDPF XB (SLU)                          840                                   0.327                        0.214
FDPF BX (SLU)                          829                                   0.346                        0.233
FDPF XB (KLU)                          864                                   0.29                         0.178
FDPF BX (KLU)                          858                                   0.305                        0.191
FDPF XB (NICSLU *)                     914                                   0.274                        0.166
FDPF BX (NICSLU *)                     903                                   0.287                        0.18
FDPF XB (CKTSO *)                      922                                   0.27                         0.163
FDPF BX (CKTSO *)                      907                                   0.286                        0.18
=====================  ======================  ===================================  ============================

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
