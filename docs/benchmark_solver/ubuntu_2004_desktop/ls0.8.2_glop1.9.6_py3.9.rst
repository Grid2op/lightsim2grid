Lightsim2grid 0.8.2 and grid2op 1.9.6 (python 3.9)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 11:59  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.9.19.final.0 (64 bit)
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
PP                                      35.1                               22                           14.4
GS                                     842                                  0.415                        0.324
GS synch                               833                                  0.425                        0.333
NR single (SLU)                       1070                                  0.161                        0.0662
NR (SLU)                              1070                                  0.162                        0.067
NR single (KLU)                       1130                                  0.111                        0.0186
NR (KLU)                              1130                                  0.11                         0.0175
NR single (NICSLU *)                  1150                                  0.109                        0.0182
NR (NICSLU *)                         1160                                  0.107                        0.0167
NR single (CKTSO *)                   1160                                  0.107                        0.017
NR (CKTSO *)                          1150                                  0.106                        0.0155
FDPF XB (SLU)                         1120                                  0.119                        0.0281
FDPF BX (SLU)                         1100                                  0.133                        0.0405
FDPF XB (KLU)                         1130                                  0.114                        0.0233
FDPF BX (KLU)                         1120                                  0.125                        0.0334
FDPF XB (NICSLU *)                    1150                                  0.113                        0.0234
FDPF BX (NICSLU *)                    1140                                  0.123                        0.0334
FDPF XB (CKTSO *)                     1140                                  0.114                        0.0231
FDPF BX (CKTSO *)                     1130                                  0.123                        0.0327
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

- date: 2024-04-23 12:05  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.9.19.final.0 (64 bit)
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
PP                                      15.9                                26.2                         17.2
GS                                       4.46                              223                          223
GS synch                                34.6                                28                           27.9
NR single (SLU)                        594                                   0.793                        0.665
NR (SLU)                               586                                   0.859                        0.651
NR single (KLU)                        980                                   0.224                        0.114
NR (KLU)                               990                                   0.213                        0.102
NR single (NICSLU *)                   976                                   0.219                        0.108
NR (NICSLU *)                          991                                   0.209                        0.0978
NR single (CKTSO *)                    991                                   0.213                        0.102
NR (CKTSO *)                           996                                   0.204                        0.0926
FDPF XB (SLU)                          905                                   0.308                        0.201
FDPF BX (SLU)                          893                                   0.324                        0.217
FDPF XB (KLU)                          936                                   0.271                        0.164
FDPF BX (KLU)                          926                                   0.285                        0.179
FDPF XB (NICSLU *)                     928                                   0.274                        0.166
FDPF BX (NICSLU *)                     907                                   0.292                        0.183
FDPF XB (CKTSO *)                      936                                   0.271                        0.164
FDPF BX (CKTSO *)                      924                                   0.285                        0.18
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
