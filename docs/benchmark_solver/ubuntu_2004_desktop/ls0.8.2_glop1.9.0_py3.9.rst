Lightsim2grid 0.8.2 and grid2op 1.9.0 (python 3.9)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 11:19  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.9.19.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.0
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
PP                                      34.8                               22                           14.5
GS                                     792                                  0.434                        0.339
GS synch                               788                                  0.441                        0.348
NR single (SLU)                       1000                                  0.164                        0.0684
NR (SLU)                              1010                                  0.166                        0.0699
NR single (KLU)                       1070                                  0.113                        0.0194
NR (KLU)                              1070                                  0.111                        0.0178
NR single (NICSLU *)                  1080                                  0.111                        0.0189
NR (NICSLU *)                         1090                                  0.108                        0.017
NR single (CKTSO *)                   1020                                  0.117                        0.0192
NR (CKTSO *)                          1030                                  0.115                        0.0174
FDPF XB (SLU)                         1060                                  0.122                        0.0293
FDPF BX (SLU)                         1040                                  0.134                        0.0421
FDPF XB (KLU)                         1060                                  0.117                        0.0243
FDPF BX (KLU)                         1050                                  0.127                        0.0347
FDPF XB (NICSLU *)                    1080                                  0.115                        0.0244
FDPF BX (NICSLU *)                    1070                                  0.126                        0.0345
FDPF XB (CKTSO *)                     1010                                  0.122                        0.0256
FDPF BX (CKTSO *)                     1000                                  0.133                        0.0367
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

- date: 2024-04-23 11:25  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.9.19.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.0
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
PP                                      16.2                                25.6                         16.9
GS                                       4.87                              204                          204
GS synch                                34.9                                27.7                         27.6
NR single (SLU)                        592                                   0.762                        0.637
NR (SLU)                               585                                   0.778                        0.654
NR single (KLU)                        859                                   0.242                        0.123
NR (KLU)                               868                                   0.228                        0.11
NR single (NICSLU *)                   860                                   0.234                        0.113
NR (NICSLU *)                          878                                   0.223                        0.104
NR single (CKTSO *)                    862                                   0.229                        0.109
NR (CKTSO *)                           867                                   0.218                        0.0988
FDPF XB (SLU)                          797                                   0.329                        0.214
FDPF BX (SLU)                          786                                   0.348                        0.233
FDPF XB (KLU)                          822                                   0.291                        0.176
FDPF BX (KLU)                          809                                   0.308                        0.193
FDPF XB (NICSLU *)                     829                                   0.288                        0.174
FDPF BX (NICSLU *)                     824                                   0.301                        0.188
FDPF XB (CKTSO *)                      818                                   0.291                        0.176
FDPF BX (CKTSO *)                      805                                   0.309                        0.194
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
