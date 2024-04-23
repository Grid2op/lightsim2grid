Lightsim2grid 0.8.2 and grid2op 1.9.2 (python 3.8)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 10:16  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
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
PP                                      31.8                               24.8                         17.2
GS                                     837                                  0.441                        0.347
GS synch                               839                                  0.444                        0.352
NR single (SLU)                       1090                                  0.166                        0.0692
NR (SLU)                              1080                                  0.168                        0.0712
NR single (KLU)                       1160                                  0.113                        0.0196
NR (KLU)                              1150                                  0.113                        0.0181
NR single (NICSLU *)                  1240                                  0.106                        0.0179
NR (NICSLU *)                         1180                                  0.112                        0.0195
NR single (CKTSO *)                   1250                                  0.103                        0.017
NR (CKTSO *)                          1250                                  0.104                        0.0155
FDPF XB (SLU)                         1140                                  0.124                        0.0307
FDPF BX (SLU)                         1130                                  0.138                        0.0438
FDPF XB (KLU)                         1150                                  0.118                        0.0251
FDPF BX (KLU)                         1140                                  0.128                        0.0361
FDPF XB (NICSLU *)                    1230                                  0.112                        0.0244
FDPF BX (NICSLU *)                    1210                                  0.122                        0.0347
FDPF XB (CKTSO *)                     1240                                  0.111                        0.0234
FDPF BX (CKTSO *)                     1230                                  0.119                        0.0336
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

- date: 2024-04-23 10:22  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
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
PP                                      15.5                                28.2                          19.5
GS                                       4.96                              201                           201
GS synch                                34.9                                27.8                          27.7
NR single (SLU)                        629                                   0.751                         0.631
NR (SLU)                               622                                   0.768                         0.648
NR single (KLU)                        928                                   0.239                         0.122
NR (KLU)                               940                                   0.226                         0.11
NR single (NICSLU *)                   989                                   0.218                         0.108
NR (NICSLU *)                          921                                   0.227                         0.108
NR single (CKTSO *)                    998                                   0.212                         0.102
NR (CKTSO *)                          1000                                   0.202                         0.092
FDPF XB (SLU)                          863                                   0.326                         0.213
FDPF BX (SLU)                          842                                   0.346                         0.233
FDPF XB (KLU)                          892                                   0.287                         0.174
FDPF BX (KLU)                          874                                   0.303                         0.192
FDPF XB (NICSLU *)                     940                                   0.271                         0.164
FDPF BX (NICSLU *)                     930                                   0.284                         0.179
FDPF XB (CKTSO *)                      943                                   0.268                         0.163
FDPF BX (CKTSO *)                      932                                   0.284                         0.178
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
