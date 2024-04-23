Lightsim2grid 0.8.2 and grid2op 1.9.1 (python 3.8)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 10:09  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.1
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
PP                                      32.1                               24.3                         16.8
GS                                     770                                  0.435                        0.342
GS synch                               769                                  0.444                        0.351
NR single (SLU)                        976                                  0.165                        0.0692
NR (SLU)                               964                                  0.167                        0.0703
NR single (KLU)                       1030                                  0.113                        0.02
NR (KLU)                              1030                                  0.113                        0.0185
NR single (NICSLU *)                  1060                                  0.11                         0.0188
NR (NICSLU *)                         1030                                  0.112                        0.0179
NR single (CKTSO *)                   1090                                  0.106                        0.0173
NR (CKTSO *)                          1080                                  0.106                        0.016
FDPF XB (SLU)                         1020                                  0.122                        0.0291
FDPF BX (SLU)                         1010                                  0.134                        0.0418
FDPF XB (KLU)                         1020                                  0.117                        0.024
FDPF BX (KLU)                         1010                                  0.127                        0.0343
FDPF XB (NICSLU *)                    1080                                  0.111                        0.023
FDPF BX (NICSLU *)                    1080                                  0.119                        0.0328
FDPF XB (CKTSO *)                     1080                                  0.111                        0.0227
FDPF BX (CKTSO *)                     1060                                  0.12                         0.0327
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

- date: 2024-04-23 10:15  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.1
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
PP                                      15                                  29                           20
GS                                       4.96                              201                          200
GS synch                                34.6                                27.9                         27.8
NR single (SLU)                        571                                   0.762                        0.635
NR (SLU)                               565                                   0.779                        0.651
NR single (KLU)                        815                                   0.244                        0.122
NR (KLU)                               825                                   0.232                        0.11
NR single (NICSLU *)                   866                                   0.224                        0.109
NR (NICSLU *)                          877                                   0.213                        0.0988
NR single (CKTSO *)                    872                                   0.219                        0.103
NR (CKTSO *)                           880                                   0.206                        0.0916
FDPF XB (SLU)                          759                                   0.333                        0.215
FDPF BX (SLU)                          747                                   0.352                        0.233
FDPF XB (KLU)                          788                                   0.291                        0.174
FDPF BX (KLU)                          775                                   0.308                        0.191
FDPF XB (NICSLU *)                     834                                   0.274                        0.164
FDPF BX (NICSLU *)                     822                                   0.291                        0.18
FDPF XB (CKTSO *)                      837                                   0.273                        0.163
FDPF BX (CKTSO *)                      818                                   0.292                        0.181
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
