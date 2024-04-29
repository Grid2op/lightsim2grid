Lightsim2grid 0.8.2 and grid2op 1.9.2 (python 3.9)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 11:32  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.9.19.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.2
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
PP                                      34.9                               22.1                         14.6
GS                                     834                                  0.435                        0.339
GS synch                               828                                  0.441                        0.344
NR single (SLU)                       1070                                  0.168                        0.0684
NR (SLU)                              1070                                  0.168                        0.069
NR single (KLU)                       1140                                  0.115                        0.0195
NR (KLU)                              1140                                  0.113                        0.0177
NR single (NICSLU *)                  1150                                  0.113                        0.0186
NR (NICSLU *)                         1150                                  0.112                        0.017
NR single (CKTSO *)                   1090                                  0.12                         0.0191
NR (CKTSO *)                          1080                                  0.119                        0.0175
FDPF XB (SLU)                         1130                                  0.124                        0.0293
FDPF BX (SLU)                         1110                                  0.137                        0.042
FDPF XB (KLU)                         1140                                  0.119                        0.0246
FDPF BX (KLU)                         1120                                  0.129                        0.0347
FDPF XB (NICSLU *)                    1120                                  0.121                        0.0252
FDPF BX (NICSLU *)                    1060                                  0.136                        0.0375
FDPF XB (CKTSO *)                     1080                                  0.126                        0.0256
FDPF BX (CKTSO *)                     1070                                  0.136                        0.0368
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

- date: 2024-04-23 11:38  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.9.19.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.2
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
PP                                       16.2                               25.8                         17.2
GS                                        5                                199                          199
GS synch                                 35                                 27.8                         27.6
NR single (SLU)                         627                                  0.754                        0.634
NR (SLU)                                620                                  0.773                        0.651
NR single (KLU)                         928                                  0.24                         0.123
NR (KLU)                                937                                  0.227                        0.11
NR single (NICSLU *)                    943                                  0.233                        0.116
NR (NICSLU *)                           947                                  0.222                        0.105
NR single (CKTSO *)                     936                                  0.227                        0.109
NR (CKTSO *)                            956                                  0.217                        0.0996
FDPF XB (SLU)                           859                                  0.328                        0.214
FDPF BX (SLU)                           843                                  0.348                        0.234
FDPF XB (KLU)                           887                                  0.29                         0.176
FDPF BX (KLU)                           873                                  0.307                        0.194
FDPF XB (NICSLU *)                      892                                  0.288                        0.175
FDPF BX (NICSLU *)                      880                                  0.304                        0.191
FDPF XB (CKTSO *)                       889                                  0.291                        0.177
FDPF BX (CKTSO *)                       875                                  0.306                        0.193
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
