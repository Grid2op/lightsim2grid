Lightsim2grid 0.8.2 and grid2op 1.9.1 (python 3.9)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 11:26  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.9.19.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.1
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
PP                                      35.8                               21.3                         14
GS                                     804                                  0.415                        0.325
GS synch                               778                                  0.434                        0.344
NR single (SLU)                        978                                  0.163                        0.0689
NR (SLU)                              1030                                  0.156                        0.066
NR single (KLU)                       1040                                  0.111                        0.0197
NR (KLU)                              1040                                  0.109                        0.0178
NR single (NICSLU *)                  1080                                  0.105                        0.0183
NR (NICSLU *)                         1100                                  0.103                        0.0164
NR single (CKTSO *)                   1110                                  0.103                        0.0171
NR (CKTSO *)                          1100                                  0.103                        0.0158
FDPF XB (SLU)                         1030                                  0.12                         0.0297
FDPF BX (SLU)                         1020                                  0.132                        0.0419
FDPF XB (KLU)                         1040                                  0.114                        0.0243
FDPF BX (KLU)                         1010                                  0.126                        0.0353
FDPF XB (NICSLU *)                    1080                                  0.109                        0.0234
FDPF BX (NICSLU *)                    1080                                  0.119                        0.0333
FDPF XB (CKTSO *)                     1090                                  0.108                        0.0228
FDPF BX (CKTSO *)                     1090                                  0.118                        0.0326
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

- date: 2024-04-23 11:31  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.9.19.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.1
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
PP                                      15.7                                26.9                         18
GS                                       4.97                              200                          200
GS synch                                34.6                                27.9                         27.7
NR single (SLU)                        573                                   0.758                        0.634
NR (SLU)                               565                                   0.777                        0.652
NR single (KLU)                        816                                   0.244                        0.123
NR (KLU)                               823                                   0.23                         0.11
NR single (NICSLU *)                   858                                   0.224                        0.108
NR (NICSLU *)                          829                                   0.223                        0.104
NR single (CKTSO *)                    876                                   0.216                        0.102
NR (CKTSO *)                           884                                   0.205                        0.0924
FDPF XB (SLU)                          761                                   0.332                        0.215
FDPF BX (SLU)                          751                                   0.351                        0.234
FDPF XB (KLU)                          786                                   0.293                        0.176
FDPF BX (KLU)                          776                                   0.31                         0.194
FDPF XB (NICSLU *)                     820                                   0.278                        0.166
FDPF BX (NICSLU *)                     818                                   0.292                        0.181
FDPF XB (CKTSO *)                      824                                   0.277                        0.167
FDPF BX (CKTSO *)                      822                                   0.29                         0.18
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
