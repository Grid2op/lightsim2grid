Lightsim2grid 0.8.2 and grid2op 1.9.8 (python 3.8)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 10:57  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.8
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
PP                                      32.5                               24.3                         16.6
GS                                     847                                  0.418                        0.328
GS synch                               816                                  0.441                        0.349
NR single (SLU)                       1050                                  0.164                        0.0683
NR (SLU)                              1080                                  0.161                        0.0671
NR single (KLU)                       1120                                  0.113                        0.0194
NR (KLU)                              1120                                  0.111                        0.0179
NR single (NICSLU *)                  1120                                  0.111                        0.0193
NR (NICSLU *)                         1140                                  0.108                        0.017
NR single (CKTSO *)                   1190                                  0.104                        0.0169
NR (CKTSO *)                          1130                                  0.109                        0.0163
FDPF XB (SLU)                         1110                                  0.121                        0.0289
FDPF BX (SLU)                         1090                                  0.133                        0.0411
FDPF XB (KLU)                         1110                                  0.117                        0.0238
FDPF BX (KLU)                         1100                                  0.127                        0.0341
FDPF XB (NICSLU *)                    1120                                  0.116                        0.0245
FDPF BX (NICSLU *)                    1100                                  0.127                        0.0351
FDPF XB (CKTSO *)                     1180                                  0.109                        0.0226
FDPF BX (CKTSO *)                     1170                                  0.119                        0.0321
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

- date: 2024-04-23 11:03  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.8
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
PP                                      14.9                                29.5                         20.3
GS                                       4.97                              200                          200
GS synch                                34.7                                27.9                         27.8
NR single (SLU)                        612                                   0.757                        0.632
NR (SLU)                               603                                   0.777                        0.649
NR single (KLU)                        890                                   0.244                        0.123
NR (KLU)                               903                                   0.23                         0.11
NR single (NICSLU *)                   913                                   0.232                        0.114
NR (NICSLU *)                          919                                   0.223                        0.104
NR single (CKTSO *)                    961                                   0.217                        0.103
NR (CKTSO *)                           973                                   0.204                        0.0915
FDPF XB (SLU)                          827                                   0.333                        0.215
FDPF BX (SLU)                          816                                   0.352                        0.236
FDPF XB (KLU)                          838                                   0.303                        0.184
FDPF BX (KLU)                          845                                   0.31                         0.193
FDPF XB (NICSLU *)                     829                                   0.305                        0.185
FDPF BX (NICSLU *)                     893                                   0.293                        0.182
FDPF XB (CKTSO *)                      912                                   0.274                        0.166
FDPF BX (CKTSO *)                      900                                   0.292                        0.181
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
