Lightsim2grid 0.8.0 and grid2op 1.10.0
======================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-03-18 14:43  CET
- system: Linux 5.15.0-56-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-6820HQ CPU @ 2.70GHz
- python version: 3.10.13.final.0 (64 bit)
- numpy version: 1.23.5
- pandas version: 2.2.1
- pandapower version: 2.13.1
- lightsim2grid version: 0.8.0
- grid2op version: 1.10.0

====================  ======================  ===================================  ============================
case14_sandbox          grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
====================  ======================  ===================================  ============================
PP                                      46.1                               18.5                          6.67
GS                                     712                                  0.539                        0.398
GS synch                               735                                  0.487                        0.345
NR single (SLU)                        916                                  0.223                        0.0809
NR (SLU)                               897                                  0.236                        0.0842
NR single (KLU)                        971                                  0.163                        0.0217
NR (KLU)                               971                                  0.162                        0.0222
NR single (NICSLU )                   979                                  0.16                         0.0212
NR (NICSLU )                          973                                  0.162                        0.0217
NR single (CKTSO )                    975                                  0.16                         0.0198
NR (CKTSO )                           964                                  0.163                        0.0212
FDPF XB (SLU)                          964                                  0.172                        0.0318
FDPF BX (SLU)                          954                                  0.183                        0.0445
FDPF XB (KLU)                          970                                  0.166                        0.0263
FDPF BX (KLU)                          964                                  0.175                        0.0371
FDPF XB (NICSLU )                     971                                  0.164                        0.026
FDPF BX (NICSLU )                     939                                  0.18                         0.0381
FDPF XB (CKTSO )                      964                                  0.169                        0.027
FDPF BX (CKTSO )                      958                                  0.177                        0.0381
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
NR single (NICSLU )                0.000122        7.63e-06          7.63e-06
NR (NICSLU )                       0.000122        7.63e-06          7.63e-06
NR single (CKTSO )                 0.000122        7.63e-06          7.63e-06
NR (CKTSO )                        0.000122        7.63e-06          7.63e-06
FDPF XB (SLU)                       0.000122        7.63e-06          7.63e-06
FDPF BX (SLU)                       0.000122        7.63e-06          7.63e-06
FDPF XB (KLU)                       0.000122        7.63e-06          7.63e-06
FDPF BX (KLU)                       0.000122        7.63e-06          7.63e-06
FDPF XB (NICSLU )                  0.000122        7.63e-06          7.63e-06
FDPF BX (NICSLU )                  0.000122        7.63e-06          7.63e-06
FDPF XB (CKTSO )                   0.000122        7.63e-06          7.63e-06
FDPF BX (CKTSO )                   0.000122        7.63e-06          7.63e-06
============================  ==============  ==============  ================


l2rpn_neurips_2020_track2_small
---------------------------------

Configuration: 

- date: 2024-03-18 14:50  CET
- system: Linux 5.15.0-56-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-6820HQ CPU @ 2.70GHz
- python version: 3.10.13.final.0 (64 bit)
- numpy version: 1.23.5
- pandas version: 2.2.1
- pandapower version: 2.13.1
- lightsim2grid version: 0.8.0
- grid2op version: 1.10.0

=====================  ======================  ===================================  ============================
neurips_2020_track2      grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
=====================  ======================  ===================================  ============================
PP                                      41.4                                20.8                           8.67
GS                                       3.56                              280                           279
GS synch                                36                                  26.7                          26.5
NR single (SLU)                        522                                   0.952                         0.776
NR (SLU)                               508                                   0.996                         0.813
NR single (KLU)                        786                                   0.314                         0.143
NR (KLU)                               787                                   0.315                         0.144
NR single (NICSLU )                   786                                   0.308                         0.135
NR (NICSLU )                          782                                   0.308                         0.136
NR single (CKTSO )                    785                                   0.301                         0.128
NR (CKTSO )                           782                                   0.306                         0.134
FDPF XB (SLU)                          736                                   0.403                         0.233
FDPF BX (SLU)                          722                                   0.424                         0.252
FDPF XB (KLU)                          753                                   0.363                         0.191
FDPF BX (KLU)                          745                                   0.381                         0.209
FDPF XB (NICSLU )                     754                                   0.363                         0.191
FDPF BX (NICSLU )                     741                                   0.38                          0.207
FDPF XB (CKTSO )                      748                                   0.364                         0.193
FDPF BX (CKTSO )                      736                                   0.382                         0.21
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
NR single (NICSLU )                      6.1e-05        0                 9.54e-07
NR (NICSLU )                             6.1e-05        0                 9.54e-07
NR single (CKTSO )                       6.1e-05        0                 9.54e-07
NR (CKTSO )                              6.1e-05        0                 9.54e-07
FDPF XB (SLU)                             6.1e-05        1.91e-06          1.53e-05
FDPF BX (SLU)                             6.1e-05        1.91e-06          7.63e-06
FDPF XB (KLU)                             6.1e-05        1.91e-06          1.53e-05
FDPF BX (KLU)                             6.1e-05        1.91e-06          7.63e-06
FDPF XB (NICSLU )                        6.1e-05        1.91e-06          1.53e-05
FDPF BX (NICSLU )                        6.1e-05        1.91e-06          7.63e-06
FDPF XB (CKTSO )                         6.1e-05        1.91e-06          1.53e-05
FDPF BX (CKTSO )                         6.1e-05        1.91e-06          7.63e-06
=================================  ==============  ==============  ================


