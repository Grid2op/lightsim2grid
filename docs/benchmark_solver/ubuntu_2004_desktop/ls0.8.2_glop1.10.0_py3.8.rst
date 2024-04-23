Lightsim2grid 0.8.2 and grid2op 1.10.0 (python 3.8)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 11:04  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
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
PP                                      36.1                               25                           17.6
GS                                     803                                  0.436                        0.342
GS synch                               793                                  0.448                        0.354
NR single (SLU)                       1020                                  0.167                        0.069
NR (SLU)                              1010                                  0.17                         0.0698
NR single (KLU)                       1080                                  0.116                        0.0197
NR (KLU)                              1080                                  0.113                        0.0182
NR single (NICSLU *)                  1080                                  0.115                        0.0194
NR (NICSLU *)                         1090                                  0.113                        0.018
NR single (CKTSO *)                   1160                                  0.106                        0.0172
NR (CKTSO *)                          1150                                  0.105                        0.0159
FDPF XB (SLU)                         1080                                  0.123                        0.029
FDPF BX (SLU)                         1060                                  0.138                        0.042
FDPF XB (KLU)                         1080                                  0.119                        0.024
FDPF BX (KLU)                         1070                                  0.129                        0.0346
FDPF XB (NICSLU *)                    1090                                  0.119                        0.0247
FDPF BX (NICSLU *)                    1090                                  0.127                        0.0346
FDPF XB (CKTSO *)                     1160                                  0.111                        0.0225
FDPF BX (CKTSO *)                     1140                                  0.121                        0.0323
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

- date: 2024-04-23 11:09  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
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
PP                                      32.4                                27.9                          20.1
GS                                       4.95                              201                           201
GS synch                                34.7                                27.9                          27.8
NR single (SLU)                        608                                   0.764                         0.643
NR (SLU)                               605                                   0.773                         0.65
NR single (KLU)                        892                                   0.24                          0.122
NR (KLU)                               909                                   0.227                         0.11
NR single (NICSLU *)                   957                                   0.22                          0.109
NR (NICSLU *)                          875                                   0.232                         0.109
NR single (CKTSO *)                    959                                   0.214                         0.103
NR (CKTSO *)                           963                                   0.203                         0.092
FDPF XB (SLU)                          833                                   0.327                         0.214
FDPF BX (SLU)                          817                                   0.348                         0.234
FDPF XB (KLU)                          862                                   0.288                         0.176
FDPF BX (KLU)                          849                                   0.304                         0.191
FDPF XB (NICSLU *)                     909                                   0.271                         0.165
FDPF BX (NICSLU *)                     904                                   0.286                         0.18
FDPF XB (CKTSO *)                      913                                   0.27                          0.163
FDPF BX (CKTSO *)                      896                                   0.286                         0.179
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
