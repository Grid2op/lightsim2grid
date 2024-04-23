Lightsim2grid 0.8.2 and grid2op 1.9.3 (python 3.8)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 10:23  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.3
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
PP                                      32.9                               23.9                         16.6
GS                                     874                                  0.419                        0.33
GS synch                               866                                  0.429                        0.34
NR single (SLU)                       1130                                  0.159                        0.0664
NR (SLU)                              1130                                  0.158                        0.0667
NR single (KLU)                       1210                                  0.108                        0.0189
NR (KLU)                              1210                                  0.108                        0.0173
NR single (NICSLU *)                  1220                                  0.106                        0.0183
NR (NICSLU *)                         1230                                  0.105                        0.0167
NR single (CKTSO *)                   1200                                  0.107                        0.0179
NR (CKTSO *)                          1200                                  0.107                        0.0165
FDPF XB (SLU)                         1190                                  0.117                        0.0279
FDPF BX (SLU)                         1170                                  0.129                        0.0399
FDPF XB (KLU)                         1200                                  0.112                        0.0232
FDPF BX (KLU)                         1190                                  0.122                        0.0329
FDPF XB (NICSLU *)                    1210                                  0.112                        0.0234
FDPF BX (NICSLU *)                    1190                                  0.121                        0.0336
FDPF XB (CKTSO *)                     1200                                  0.112                        0.0233
FDPF BX (CKTSO *)                     1140                                  0.127                        0.0347
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

- date: 2024-04-23 10:29  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.3
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
PP                                      15                                  29.5                         20.5
GS                                       4.95                              201                          201
GS synch                                34.8                                27.8                         27.7
NR single (SLU)                        620                                   0.756                        0.632
NR (SLU)                               613                                   0.774                        0.65
NR single (KLU)                        910                                   0.244                        0.123
NR (KLU)                               920                                   0.231                        0.112
NR single (NICSLU *)                   931                                   0.232                        0.114
NR (NICSLU *)                          938                                   0.222                        0.104
NR single (CKTSO *)                    975                                   0.217                        0.104
NR (CKTSO *)                           986                                   0.205                        0.0918
FDPF XB (SLU)                          847                                   0.328                        0.213
FDPF BX (SLU)                          830                                   0.349                        0.233
FDPF XB (KLU)                          878                                   0.289                        0.174
FDPF BX (KLU)                          860                                   0.305                        0.19
FDPF XB (NICSLU *)                     927                                   0.274                        0.164
FDPF BX (NICSLU *)                     910                                   0.29                         0.181
FDPF XB (CKTSO *)                      925                                   0.273                        0.164
FDPF BX (CKTSO *)                      915                                   0.287                        0.179
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
