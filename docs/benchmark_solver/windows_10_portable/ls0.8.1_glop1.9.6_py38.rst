Lightsim2grid 0.8.1 and grid2op 1.9.6 (python 3.8)
====================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-03-28 07:21  Paris, Madrid
- system: Windows 10
- OS:
- processor: AMD Ryzen 7 4800HS with Radeon Graphics
- python version: 3.8.3.final.0 (64 bit)
- numpy version: 1.21.0
- pandas version: 2.0.3
- pandapower version: 2.13.1
- grid2op version: 1.9.6
- lightsim2grid version: 0.8.1
- lightsim2grid extra information:

        - klu_solver_available: True
        - nicslu_solver_available: False
        - cktso_solver_available: False
        - compiled_march_native: True
        - compiled_o3_optim: True

================  ======================  ===================================  ============================
case14_sandbox      grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
================  ======================  ===================================  ============================
PP                                    49                               14.6                          5.09
GS                                   733                                0.593                        0.495
GS synch                             653                                0.716                        0.611
NR single (SLU)                     1060                                0.168                        0.0674
NR (SLU)                            1050                                0.172                        0.0715
NR single (KLU)                     1130                                0.118                        0.0205
NR (KLU)                            1130                                0.12                         0.0233
FDPF XB (SLU)                       1110                                0.136                        0.0393
FDPF BX (SLU)                       1080                                0.155                        0.0573
FDPF XB (KLU)                       1120                                0.126                        0.0288
FDPF BX (KLU)                       1090                                0.14                         0.0418
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

- date: 2024-03-28 07:28  Paris, Madrid
- system: Windows 10
- OS:
- processor: AMD Ryzen 7 4800HS with Radeon Graphics
- python version: 3.8.3.final.0 (64 bit)
- numpy version: 1.21.0
- pandas version: 2.0.3
- pandapower version: 2.13.1
- grid2op version: 1.9.6
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
PP                                      21.3                                16.9                           6.58
GS                                       6.26                              158                           158
GS synch                                26.7                                36.2                          36
NR single (SLU)                        556                                   0.767                         0.609
NR (SLU)                               592                                   0.769                         0.639
NR single (KLU)                        803                                   0.256                         0.114
NR (KLU)                               849                                   0.257                         0.128
FDPF XB (SLU)                          717                                   0.416                         0.278
FDPF BX (SLU)                          735                                   0.43                          0.301
FDPF XB (KLU)                          708                                   0.367                         0.215
FDPF BX (KLU)                          734                                   0.373                         0.234
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
