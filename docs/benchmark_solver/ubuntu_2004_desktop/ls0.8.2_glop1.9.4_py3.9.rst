Lightsim2grid 0.8.2 and grid2op 1.9.4 (python 3.9)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 11:45  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.9.19.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.4
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
PP                                      35.2                               21.8                         14.3
GS                                     854                                  0.43                         0.338
GS synch                               849                                  0.438                        0.347
NR single (SLU)                       1090                                  0.168                        0.0702
NR (SLU)                              1090                                  0.169                        0.0705
NR single (KLU)                       1170                                  0.113                        0.0195
NR (KLU)                              1160                                  0.112                        0.0186
NR single (NICSLU *)                  1210                                  0.109                        0.0187
NR (NICSLU *)                         1200                                  0.108                        0.0171
NR single (CKTSO *)                   1120                                  0.116                        0.0193
NR (CKTSO *)                          1130                                  0.114                        0.0175
FDPF XB (SLU)                         1160                                  0.124                        0.0297
FDPF BX (SLU)                         1140                                  0.135                        0.0416
FDPF XB (KLU)                         1170                                  0.116                        0.0242
FDPF BX (KLU)                         1150                                  0.128                        0.035
FDPF XB (NICSLU *)                    1170                                  0.116                        0.0248
FDPF BX (NICSLU *)                    1100                                  0.135                        0.0373
FDPF XB (CKTSO *)                     1120                                  0.121                        0.0253
FDPF BX (CKTSO *)                     1110                                  0.132                        0.0366
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

- date: 2024-04-23 11:51  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.9.19.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.4
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
PP                                      15.5                                26.5                         17.5
GS                                       4.78                              208                          208
GS synch                                34.8                                27.9                         27.8
NR single (SLU)                        627                                   0.76                         0.636
NR (SLU)                               622                                   0.778                        0.654
NR single (KLU)                        937                                   0.242                        0.124
NR (KLU)                               950                                   0.228                        0.11
NR single (NICSLU *)                  1000                                   0.22                         0.108
NR (NICSLU *)                         1010                                   0.211                        0.0991
NR single (CKTSO *)                   1020                                   0.211                        0.102
NR (CKTSO *)                          1030                                   0.203                        0.0922
FDPF XB (SLU)                          865                                   0.331                        0.217
FDPF BX (SLU)                          852                                   0.351                        0.236
FDPF XB (KLU)                          897                                   0.292                        0.178
FDPF BX (KLU)                          886                                   0.307                        0.194
FDPF XB (NICSLU *)                     953                                   0.273                        0.166
FDPF BX (NICSLU *)                     938                                   0.289                        0.181
FDPF XB (CKTSO *)                      886                                   0.357                        0.165
FDPF BX (CKTSO *)                      944                                   0.288                        0.18
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
