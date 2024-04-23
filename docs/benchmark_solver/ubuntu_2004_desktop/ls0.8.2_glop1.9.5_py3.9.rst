Lightsim2grid 0.8.2 and grid2op 1.9.5 (python 3.9)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 11:52  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.9.19.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.5
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
PP                                      35.1                               22                           14.5
GS                                     873                                  0.416                        0.325
GS synch                               862                                  0.428                        0.336
NR single (SLU)                       1110                                  0.162                        0.0674
NR (SLU)                              1120                                  0.161                        0.0673
NR single (KLU)                       1190                                  0.11                         0.0191
NR (KLU)                              1190                                  0.109                        0.0174
NR single (NICSLU *)                  1160                                  0.113                        0.0193
NR (NICSLU *)                         1190                                  0.109                        0.0172
NR single (CKTSO *)                   1120                                  0.117                        0.0193
NR (CKTSO *)                          1110                                  0.115                        0.0175
FDPF XB (SLU)                         1180                                  0.118                        0.0284
FDPF BX (SLU)                         1160                                  0.13                         0.0404
FDPF XB (KLU)                         1190                                  0.113                        0.0237
FDPF BX (KLU)                         1180                                  0.123                        0.0337
FDPF XB (NICSLU *)                    1160                                  0.116                        0.0252
FDPF BX (NICSLU *)                    1150                                  0.126                        0.0351
FDPF XB (CKTSO *)                     1100                                  0.123                        0.026
FDPF BX (CKTSO *)                     1100                                  0.132                        0.0366
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

- date: 2024-04-23 11:58  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.9.19.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.5
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
PP                                      16.2                                25.7                         17.1
GS                                       4.98                              200                          200
GS synch                                35                                  27.7                         27.6
NR single (SLU)                        632                                   0.756                        0.634
NR (SLU)                               624                                   0.775                        0.652
NR single (KLU)                        944                                   0.239                        0.122
NR (KLU)                               956                                   0.226                        0.109
NR single (NICSLU *)                   952                                   0.23                         0.114
NR (NICSLU *)                          958                                   0.222                        0.105
NR single (CKTSO *)                   1010                                   0.215                        0.104
NR (CKTSO *)                          1020                                   0.203                        0.0923
FDPF XB (SLU)                          873                                   0.326                        0.213
FDPF BX (SLU)                          857                                   0.347                        0.233
FDPF XB (KLU)                          900                                   0.288                        0.176
FDPF BX (KLU)                          889                                   0.305                        0.192
FDPF XB (NICSLU *)                     895                                   0.291                        0.177
FDPF BX (NICSLU *)                     882                                   0.307                        0.193
FDPF XB (CKTSO *)                      875                                   0.359                        0.167
FDPF BX (CKTSO *)                      935                                   0.289                        0.182
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
