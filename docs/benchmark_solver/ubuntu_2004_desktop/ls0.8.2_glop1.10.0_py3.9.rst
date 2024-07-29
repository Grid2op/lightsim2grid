Lightsim2grid 0.8.2 and grid2op 1.10.0 (python 3.9)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 12:19  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.9.19.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.10.0
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
PP                                      42.2                               21                           14
GS                                     844                                  0.424                        0.333
GS synch                               842                                  0.426                        0.335
NR single (SLU)                       1080                                  0.16                         0.0661
NR (SLU)                              1030                                  0.208                        0.0675
NR single (KLU)                       1150                                  0.11                         0.019
NR (KLU)                              1150                                  0.108                        0.0174
NR single (NICSLU *)                  1140                                  0.112                        0.0188
NR (NICSLU *)                         1140                                  0.109                        0.0178
NR single (CKTSO *)                   1100                                  0.113                        0.0185
NR (CKTSO *)                          1090                                  0.114                        0.0174
FDPF XB (SLU)                         1130                                  0.12                         0.0288
FDPF BX (SLU)                         1120                                  0.131                        0.0409
FDPF XB (KLU)                         1140                                  0.114                        0.0237
FDPF BX (KLU)                         1130                                  0.124                        0.0338
FDPF XB (NICSLU *)                    1140                                  0.115                        0.0244
FDPF BX (NICSLU *)                    1080                                  0.131                        0.0365
FDPF XB (CKTSO *)                     1090                                  0.119                        0.0248
FDPF BX (CKTSO *)                     1080                                  0.13                         0.0356
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

- date: 2024-04-23 12:25  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.9.19.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.10.0
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
PP                                      36.6                                24.4                         17
GS                                       4.59                              217                          217
GS synch                                34.6                                28                           27.8
NR single (SLU)                        577                                   0.799                        0.666
NR (SLU)                               601                                   0.776                        0.65
NR single (KLU)                        944                                   0.229                        0.115
NR (KLU)                               953                                   0.218                        0.103
NR single (NICSLU *)                   951                                   0.221                        0.107
NR (NICSLU *)                          958                                   0.213                        0.0976
NR single (CKTSO *)                    953                                   0.216                        0.102
NR (CKTSO *)                           956                                   0.208                        0.0929
FDPF XB (SLU)                          878                                   0.31                         0.199
FDPF BX (SLU)                          863                                   0.328                        0.217
FDPF XB (KLU)                          906                                   0.276                        0.165
FDPF BX (KLU)                          892                                   0.291                        0.18
FDPF XB (NICSLU *)                     898                                   0.278                        0.166
FDPF BX (NICSLU *)                     884                                   0.294                        0.182
FDPF XB (CKTSO *)                      908                                   0.275                        0.164
FDPF BX (CKTSO *)                      889                                   0.293                        0.181
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
