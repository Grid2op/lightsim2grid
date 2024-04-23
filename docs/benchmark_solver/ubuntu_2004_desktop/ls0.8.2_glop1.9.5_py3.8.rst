Lightsim2grid 0.8.2 and grid2op 1.9.5 (python 3.8)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 10:36  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.5
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
PP                                      32.6                               24.1                         16.6
GS                                     862                                  0.42                         0.331
GS synch                               853                                  0.432                        0.341
NR single (SLU)                       1070                                  0.166                        0.069
NR (SLU)                              1100                                  0.162                        0.0673
NR single (KLU)                       1150                                  0.112                        0.0194
NR (KLU)                              1150                                  0.111                        0.018
NR single (NICSLU *)                  1220                                  0.107                        0.0181
NR (NICSLU *)                         1220                                  0.104                        0.0166
NR single (CKTSO *)                   1230                                  0.104                        0.017
NR (CKTSO *)                          1220                                  0.103                        0.0156
FDPF XB (SLU)                         1130                                  0.123                        0.0293
FDPF BX (SLU)                         1120                                  0.135                        0.0414
FDPF XB (KLU)                         1150                                  0.114                        0.0237
FDPF BX (KLU)                         1130                                  0.126                        0.0344
FDPF XB (NICSLU *)                    1200                                  0.11                         0.023
FDPF BX (NICSLU *)                    1200                                  0.119                        0.0328
FDPF XB (CKTSO *)                     1220                                  0.109                        0.0223
FDPF BX (CKTSO *)                     1210                                  0.117                        0.032
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

- date: 2024-04-23 10:42  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.5
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
PP                                      15.2                                28.5                         19.6
GS                                       4.97                              200                          200
GS synch                                34.8                                27.9                         27.8
NR single (SLU)                        635                                   0.754                        0.636
NR (SLU)                               626                                   0.77                         0.651
NR single (KLU)                        939                                   0.242                        0.127
NR (KLU)                               947                                   0.229                        0.112
NR single (NICSLU *)                   967                                   0.229                        0.118
NR (NICSLU *)                         1010                                   0.209                        0.0994
NR single (CKTSO *)                   1020                                   0.214                        0.106
NR (CKTSO *)                          1020                                   0.201                        0.0924
FDPF XB (SLU)                          873                                   0.325                        0.214
FDPF BX (SLU)                          861                                   0.344                        0.234
FDPF XB (KLU)                          905                                   0.285                        0.175
FDPF BX (KLU)                          896                                   0.3                          0.191
FDPF XB (NICSLU *)                     961                                   0.268                        0.164
FDPF BX (NICSLU *)                     946                                   0.283                        0.179
FDPF XB (CKTSO *)                      953                                   0.271                        0.167
FDPF BX (CKTSO *)                      948                                   0.282                        0.179
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
