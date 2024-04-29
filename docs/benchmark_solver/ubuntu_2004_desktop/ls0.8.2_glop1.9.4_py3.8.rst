Lightsim2grid 0.8.2 and grid2op 1.9.4 (python 3.8)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 10:29  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
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
PP                                      32.1                               24.4                         16.9
GS                                     851                                  0.435                        0.342
GS synch                               841                                  0.445                        0.352
NR single (SLU)                       1090                                  0.166                        0.0695
NR (SLU)                              1090                                  0.167                        0.0701
NR single (KLU)                       1170                                  0.114                        0.0198
NR (KLU)                              1160                                  0.113                        0.0188
NR single (NICSLU *)                  1190                                  0.111                        0.0193
NR (NICSLU *)                         1200                                  0.109                        0.0174
NR single (CKTSO *)                   1170                                  0.112                        0.0184
NR (CKTSO *)                          1170                                  0.111                        0.0169
FDPF XB (SLU)                         1160                                  0.122                        0.0294
FDPF BX (SLU)                         1140                                  0.134                        0.0418
FDPF XB (KLU)                         1170                                  0.117                        0.0239
FDPF BX (KLU)                         1150                                  0.129                        0.0348
FDPF XB (NICSLU *)                    1160                                  0.118                        0.0246
FDPF BX (NICSLU *)                    1170                                  0.125                        0.0346
FDPF XB (CKTSO *)                     1160                                  0.118                        0.0242
FDPF BX (CKTSO *)                     1150                                  0.127                        0.0346
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

- date: 2024-04-23 10:35  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
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
PP                                      14.9                                29.7                         20.4
GS                                       4.68                              213                          213
GS synch                                34.7                                28                           27.8
NR single (SLU)                        625                                   0.757                        0.636
NR (SLU)                               574                                   0.813                        0.68
NR single (KLU)                        919                                   0.242                        0.125
NR (KLU)                               933                                   0.228                        0.111
NR single (NICSLU *)                   983                                   0.218                        0.108
NR (NICSLU *)                          995                                   0.208                        0.0977
NR single (CKTSO *)                    995                                   0.212                        0.103
NR (CKTSO *)                          1010                                   0.202                        0.0919
FDPF XB (SLU)                          858                                   0.326                        0.213
FDPF BX (SLU)                          842                                   0.346                        0.232
FDPF XB (KLU)                          887                                   0.286                        0.174
FDPF BX (KLU)                          872                                   0.303                        0.19
FDPF XB (NICSLU *)                     932                                   0.273                        0.166
FDPF BX (NICSLU *)                     924                                   0.285                        0.179
FDPF XB (CKTSO *)                      942                                   0.27                         0.164
FDPF BX (CKTSO *)                      927                                   0.285                        0.179
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
