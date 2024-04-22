Lightsim2grid 0.8.0 and grid2op 1.9.8
======================================

l2rpn_case14_sandbox
----------------------

Configuration:

- date: 2024-03-18 14:52  CET
- system: Linux 5.15.0-56-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-6820HQ CPU @ 2.70GHz
- python version: 3.10.13.final.0 (64 bit)
- numpy version: 1.23.5
- pandas version: 2.2.1
- pandapower version: 2.13.1
- lightsim2grid version: 0.8.0
- grid2op version: 1.9.8

====================  ======================  ===================================  ============================
case14_sandbox          grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
====================  ======================  ===================================  ============================
PP                                      37.9                               18.7                          6.43
GS                                     728                                  0.535                        0.395
GS synch                               758                                  0.483                        0.345
NR single (SLU)                        952                                  0.219                        0.0798
NR (SLU)                               929                                  0.231                        0.0843
NR single (KLU)                        979                                  0.164                        0.0222
NR (KLU)                              1010                                  0.159                        0.022
NR single (NICSLU)                  1000                                  0.16                         0.0213
NR (NICSLU)                         1010                                  0.159                        0.0218
NR single (CKTSO)                    992                                  0.159                        0.0204
NR (CKTSO)                          1000                                  0.159                        0.0212
FDPF XB (SLU)                         1010                                  0.169                        0.0323
FDPF BX (SLU)                          989                                  0.182                        0.0455
FDPF XB (KLU)                         1010                                  0.164                        0.0263
FDPF BX (KLU)                         1000                                  0.174                        0.0376
FDPF XB (NICSLU)                     999                                  0.164                        0.0263
FDPF BX (NICSLU)                     985                                  0.176                        0.0379
FDPF XB (CKTSO)                      994                                  0.166                        0.0266
FDPF BX (CKTSO)                      990                                  0.173                        0.0374
====================  ======================  ===================================  ============================


Differences:

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
NR single (NICSLU)                0.000122        7.63e-06          7.63e-06
NR (NICSLU)                       0.000122        7.63e-06          7.63e-06
NR single (CKTSO)                 0.000122        7.63e-06          7.63e-06
NR (CKTSO)                        0.000122        7.63e-06          7.63e-06
FDPF XB (SLU)                       0.000122        7.63e-06          7.63e-06
FDPF BX (SLU)                       0.000122        7.63e-06          7.63e-06
FDPF XB (KLU)                       0.000122        7.63e-06          7.63e-06
FDPF BX (KLU)                       0.000122        7.63e-06          7.63e-06
FDPF XB (NICSLU)                  0.000122        7.63e-06          7.63e-06
FDPF BX (NICSLU)                  0.000122        7.63e-06          7.63e-06
FDPF XB (CKTSO)                   0.000122        7.63e-06          7.63e-06
FDPF BX (CKTSO)                   0.000122        7.63e-06          7.63e-06
============================  ==============  ==============  ================

l2rpn_neurips_2020_track2_small
--------------------------------

Configuration: 

- date: 2024-03-18 15:00  CET
- system: Linux 5.15.0-56-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-6820HQ CPU @ 2.70GHz
- python version: 3.10.13.final.0 (64 bit)
- numpy version: 1.23.5
- pandas version: 2.2.1
- pandapower version: 2.13.1
- lightsim2grid version: 0.8.0
- grid2op version: 1.9.8

=====================  ======================  ===================================  ============================
neurips_2020_track2      grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
=====================  ======================  ===================================  ============================
PP                                       16.3                               21.2                           8.24
GS                                        3.5                              284                           284
GS synch                                 35.7                               26.9                          26.7
NR single (SLU)                         531                                  0.952                         0.776
NR (SLU)                                513                                  0.995                         0.812
NR single (KLU)                         798                                  0.317                         0.145
NR (KLU)                                796                                  0.317                         0.145
NR single (NICSLU)                    806                                  0.307                         0.136
NR (NICSLU)                           812                                  0.305                         0.134
NR single (CKTSO)                     815                                  0.299                         0.129
NR (CKTSO)                            809                                  0.302                         0.131
FDPF XB (SLU)                           736                                  0.415                         0.24
FDPF BX (SLU)                           735                                  0.425                         0.254
FDPF XB (KLU)                           769                                  0.366                         0.195
FDPF BX (KLU)                           761                                  0.381                         0.211
FDPF XB (NICSLU)                      774                                  0.362                         0.191
FDPF BX (NICSLU)                      766                                  0.378                         0.209
FDPF XB (CKTSO)                       736                                  0.38                          0.202
FDPF BX (CKTSO)                       729                                  0.4                           0.222
=====================  ======================  ===================================  ============================


Differences:
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
NR single (NICSLU)                      6.1e-05        0                 9.54e-07
NR (NICSLU)                             6.1e-05        0                 9.54e-07
NR single (CKTSO)                       6.1e-05        0                 9.54e-07
NR (CKTSO)                              6.1e-05        0                 9.54e-07
FDPF XB (SLU)                             6.1e-05        1.91e-06          1.53e-05
FDPF BX (SLU)                             6.1e-05        1.91e-06          7.63e-06
FDPF XB (KLU)                             6.1e-05        1.91e-06          1.53e-05
FDPF BX (KLU)                             6.1e-05        1.91e-06          7.63e-06
FDPF XB (NICSLU)                        6.1e-05        1.91e-06          1.53e-05
FDPF BX (NICSLU)                        6.1e-05        1.91e-06          7.63e-06
FDPF XB (CKTSO)                         6.1e-05        1.91e-06          1.53e-05
FDPF BX (CKTSO)                         6.1e-05        1.91e-06          7.63e-06
=================================  ==============  ==============  ================

