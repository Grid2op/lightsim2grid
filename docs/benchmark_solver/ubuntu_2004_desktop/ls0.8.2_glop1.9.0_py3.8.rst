Lightsim2grid 0.8.2 and grid2op 1.9.0 (python 3.8)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 10:02  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.0
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
PP                                      31.7                               24.7                         17.1
GS                                     796                                  0.432                        0.34
GS synch                               784                                  0.449                        0.357
NR single (SLU)                       1000                                  0.166                        0.0688
NR (SLU)                              1000                                  0.167                        0.0698
NR single (KLU)                       1070                                  0.111                        0.0198
NR (KLU)                              1070                                  0.11                         0.0184
NR single (NICSLU *)                  1060                                  0.111                        0.0194
NR (NICSLU *)                         1090                                  0.109                        0.0174
NR single (CKTSO *)                   1150                                  0.102                        0.017
NR (CKTSO *)                          1140                                  0.102                        0.016
FDPF XB (SLU)                         1060                                  0.121                        0.0291
FDPF BX (SLU)                         1040                                  0.134                        0.0415
FDPF XB (KLU)                         1070                                  0.115                        0.0243
FDPF BX (KLU)                         1060                                  0.125                        0.0347
FDPF XB (NICSLU *)                    1010                                  0.121                        0.0256
FDPF BX (NICSLU *)                    1000                                  0.132                        0.037
FDPF XB (CKTSO *)                     1140                                  0.107                        0.0223
FDPF BX (CKTSO *)                     1120                                  0.118                        0.032
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

- date: 2024-04-23 10:08  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.0
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
PP                                      15.2                                28.7                         19.7
GS                                       4.09                              243                          243
GS synch                                34.2                                28.3                         28.2
NR single (SLU)                        632                                   0.704                        0.59
NR (SLU)                               587                                   0.772                        0.649
NR single (KLU)                        914                                   0.226                        0.115
NR (KLU)                               919                                   0.214                        0.103
NR single (NICSLU *)                   876                                   0.228                        0.113
NR (NICSLU *)                          922                                   0.209                        0.0979
NR single (CKTSO *)                    930                                   0.213                        0.102
NR (CKTSO *)                           931                                   0.203                        0.0923
FDPF XB (SLU)                          847                                   0.308                        0.201
FDPF BX (SLU)                          836                                   0.325                        0.218
FDPF XB (KLU)                          875                                   0.27                         0.163
FDPF BX (KLU)                          868                                   0.284                        0.178
FDPF XB (NICSLU *)                     822                                   0.29                         0.176
FDPF BX (NICSLU *)                     859                                   0.288                        0.181
FDPF XB (CKTSO *)                      882                                   0.269                        0.163
FDPF BX (CKTSO *)                      865                                   0.287                        0.18
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
