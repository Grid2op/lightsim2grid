Lightsim2grid 0.8.2 and grid2op 1.9.8 (python 3.9)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 12:13  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.9.19.final.0 (64 bit)
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
PP                                      34.9                               22.1                         14.5
GS                                     852                                  0.42                         0.327
GS synch                               844                                  0.43                         0.338
NR single (SLU)                       1090                                  0.162                        0.0668
NR (SLU)                              1090                                  0.163                        0.0674
NR single (KLU)                       1150                                  0.112                        0.0188
NR (KLU)                              1150                                  0.111                        0.0177
NR single (NICSLU *)                  1160                                  0.111                        0.0184
NR (NICSLU *)                         1180                                  0.108                        0.0164
NR single (CKTSO *)                   1160                                  0.111                        0.0178
NR (CKTSO *)                          1160                                  0.11                         0.0161
FDPF XB (SLU)                         1150                                  0.121                        0.0286
FDPF BX (SLU)                         1120                                  0.134                        0.0407
FDPF XB (KLU)                         1150                                  0.116                        0.0236
FDPF BX (KLU)                         1130                                  0.127                        0.0338
FDPF XB (NICSLU *)                    1160                                  0.116                        0.0239
FDPF BX (NICSLU *)                    1150                                  0.125                        0.0341
FDPF XB (CKTSO *)                     1150                                  0.116                        0.0238
FDPF BX (CKTSO *)                     1130                                  0.127                        0.0341
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

- date: 2024-04-23 12:18  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.9.19.final.0 (64 bit)
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
PP                                      15.8                                26.1                         17.1
GS                                       4.93                              202                          202
GS synch                                34.3                                28.3                         28.2
NR single (SLU)                        620                                   0.756                        0.634
NR (SLU)                               582                                   0.859                        0.651
NR single (KLU)                        909                                   0.242                        0.124
NR (KLU)                               920                                   0.229                        0.111
NR single (NICSLU *)                   964                                   0.222                        0.109
NR (NICSLU *)                          983                                   0.21                         0.098
NR single (CKTSO *)                    979                                   0.215                        0.102
NR (CKTSO *)                           959                                   0.212                        0.0978
FDPF XB (SLU)                          841                                   0.332                        0.217
FDPF BX (SLU)                          830                                   0.351                        0.236
FDPF XB (KLU)                          873                                   0.292                        0.179
FDPF BX (KLU)                          853                                   0.311                        0.197
FDPF XB (NICSLU *)                     928                                   0.274                        0.167
FDPF BX (NICSLU *)                     908                                   0.293                        0.183
FDPF XB (CKTSO *)                      917                                   0.278                        0.17
FDPF BX (CKTSO *)                      910                                   0.29                         0.182
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
