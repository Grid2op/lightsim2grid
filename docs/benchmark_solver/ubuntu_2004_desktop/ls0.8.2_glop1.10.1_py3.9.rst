Lightsim2grid 0.8.2 and grid2op 1.10.1 (python 3.9)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 12:25  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.9.19.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.10.1
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
PP                                      40.3                               22                           14.8
GS                                     803                                  0.434                        0.338
GS synch                               799                                  0.44                         0.343
NR single (SLU)                       1020                                  0.169                        0.0688
NR (SLU)                               982                                  0.216                        0.0683
NR single (KLU)                       1080                                  0.117                        0.0194
NR (KLU)                              1080                                  0.115                        0.0179
NR single (NICSLU *)                  1030                                  0.123                        0.02
NR (NICSLU *)                         1080                                  0.114                        0.0175
NR single (CKTSO *)                   1160                                  0.108                        0.0167
NR (CKTSO *)                          1150                                  0.108                        0.0156
FDPF XB (SLU)                         1070                                  0.126                        0.029
FDPF BX (SLU)                         1060                                  0.138                        0.0413
FDPF XB (KLU)                         1080                                  0.12                         0.024
FDPF BX (KLU)                         1060                                  0.132                        0.0347
FDPF XB (NICSLU *)                    1140                                  0.114                        0.0233
FDPF BX (NICSLU *)                    1130                                  0.123                        0.0332
FDPF XB (CKTSO *)                     1150                                  0.113                        0.0227
FDPF BX (CKTSO *)                     1130                                  0.124                        0.0324
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

- date: 2024-04-23 12:31  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.9.19.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.10.1
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
PP                                       34.2                               26.1                         18.2
GS                                        4.6                              216                          216
GS synch                                 34.8                               27.8                         27.7
NR single (SLU)                         600                                  0.758                        0.633
NR (SLU)                                594                                  0.774                        0.649
NR single (KLU)                         871                                  0.244                        0.122
NR (KLU)                                881                                  0.23                         0.109
NR single (NICSLU *)                    928                                  0.223                        0.108
NR (NICSLU *)                           882                                  0.227                        0.105
NR single (CKTSO *)                     932                                  0.217                        0.102
NR (CKTSO *)                            941                                  0.207                        0.0932
FDPF XB (SLU)                           809                                  0.33                         0.213
FDPF BX (SLU)                           793                                  0.352                        0.233
FDPF XB (KLU)                           835                                  0.292                        0.176
FDPF BX (KLU)                           824                                  0.309                        0.192
FDPF XB (NICSLU *)                      883                                  0.276                        0.166
FDPF BX (NICSLU *)                      875                                  0.29                         0.181
FDPF XB (CKTSO *)                       886                                  0.275                        0.165
FDPF BX (CKTSO *)                       875                                  0.29                         0.18
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
