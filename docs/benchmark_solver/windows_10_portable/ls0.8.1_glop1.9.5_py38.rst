Lightsim2grid 0.8.1 and grid2op 1.9.5 (python 3.8)
====================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-03-27 22:19  Paris, Madrid
- system: Windows 10
- OS:
- processor: AMD Ryzen 7 4800HS with Radeon Graphics
- python version: 3.8.3.final.0 (64 bit)
- numpy version: 1.21.0
- pandas version: 2.0.3
- pandapower version: 2.13.1
- grid2op version: 1.9.5
- lightsim2grid version: 0.8.1
- lightsim2grid extra information:

        - klu_solver_available: True
        - nicslu_solver_available: False
        - cktso_solver_available: False
        - compiled_march_native: True
        - compiled_o3_optim: True

**NB** Due to issue with dll, cktso and nicslu could not be benchmarked

================  ======================  ===================================  ============================
case14_sandbox      grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
================  ======================  ===================================  ============================
PP                                  48.1                               14.9                          5.16
GS                                 728                                  0.6                          0.501
GS synch                           663                                  0.729                        0.63
NR single (SLU)                   1070                                  0.167                        0.0673
NR (SLU)                          1050                                  0.172                        0.0714
NR single (KLU)                   1120                                  0.119                        0.0212
NR (KLU)                          1120                                  0.122                        0.0237
FDPF XB (SLU)                     1100                                  0.136                        0.0393
FDPF BX (SLU)                     1070                                  0.156                        0.0572
FDPF XB (KLU)                     1110                                  0.128                        0.0293
FDPF BX (KLU)                     1090                                  0.14                         0.0424
================  ======================  ===================================  ============================

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
FDPF XB (SLU)                       0.000122        7.63e-06          7.63e-06
FDPF BX (SLU)                       0.000122        7.63e-06          7.63e-06
FDPF XB (KLU)                       0.000122        7.63e-06          7.63e-06
FDPF BX (KLU)                       0.000122        7.63e-06          7.63e-06
============================  ==============  ==============  ================

l2rpn_neurips_2020_track2_small
---------------------------------

Configuration:

- date: 2024-03-27 22:28  Paris, Madrid
- system: Windows 10
- OS:
- processor: AMD Ryzen 7 4800HS with Radeon Graphics
- python version: 3.8.3.final.0 (64 bit)
- numpy version: 1.21.0
- pandas version: 2.0.3
- pandapower version: 2.13.1
- grid2op version: 1.9.5
- lightsim2grid version: 0.8.1
- lightsim2grid extra information:

        - klu_solver_available: True
        - nicslu_solver_available: False
        - cktso_solver_available: False
        - compiled_march_native: True
        - compiled_o3_optim: True

=====================  ======================  ===================================  ============================
neurips_2020_track2      grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
=====================  ======================  ===================================  ============================
PP                                      21.2                                17                             6.61
GS                                       6.28                              158                           158
GS synch                                26.8                                36.1                          35.9
NR single (SLU)                        616                                   0.716                         0.583
NR (SLU)                               564                                   0.789                         0.637
NR single (KLU)                        921                                   0.228                         0.105
NR (KLU)                               901                                   0.245                         0.122
FDPF XB (SLU)                          795                                   0.389                         0.268
FDPF BX (SLU)                          783                                   0.414                         0.293
FDPF XB (KLU)                          848                                   0.323                         0.203
FDPF BX (KLU)                          838                                   0.34                          0.222
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
FDPF XB (SLU)                             6.1e-05        1.91e-06          1.53e-05
FDPF BX (SLU)                             6.1e-05        1.91e-06          7.63e-06
FDPF XB (KLU)                             6.1e-05        1.91e-06          1.53e-05
FDPF BX (KLU)                             6.1e-05        1.91e-06          7.63e-06
=================================  ==============  ==============  ================