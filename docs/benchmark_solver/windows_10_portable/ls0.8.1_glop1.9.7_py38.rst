Lightsim2grid 0.8.1 and grid2op 1.9.7 (python 3.8)
====================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-03-28 07:56  Paris, Madrid
- system: Windows 10
- OS:
- processor: AMD Ryzen 7 4800HS with Radeon Graphics
- python version: 3.8.3.final.0 (64 bit)
- numpy version: 1.21.0
- pandas version: 2.0.3
- pandapower version: 2.13.1
- grid2op version: 1.9.7
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
PP                                  49.4                               14.5                          5.04
GS                                 734                                  0.591                        0.492
GS synch                           676                                  0.704                        0.605
NR single (SLU)                   1060                                  0.166                        0.0663
NR (SLU)                          1050                                  0.172                        0.0715
NR single (KLU)                   1130                                  0.118                        0.0202
NR (KLU)                          1120                                  0.122                        0.023
FDPF XB (SLU)                     1090                                  0.138                        0.0388
FDPF BX (SLU)                     1090                                  0.153                        0.0556
FDPF XB (KLU)                     1120                                  0.126                        0.0283
FDPF BX (KLU)                     1100                                  0.138                        0.0411
================  ======================  ===================================  ============================

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

- date: 2024-03-28 08:15  Paris, Madrid
- system: Windows 10
- OS:
- processor: AMD Ryzen 7 4800HS with Radeon Graphics
- python version: 3.8.3.final.0 (64 bit)
- numpy version: 1.21.0
- pandas version: 2.0.3
- pandapower version: 2.13.1
- grid2op version: 1.9.7
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
PP                                      21.4                                16.9                           6.57
GS                                       6.31                              157                           157
GS synch                                26.7                                36.3                          36.1
NR single (SLU)                        607                                   0.729                         0.595
NR (SLU)                               593                                   0.769                         0.635
NR single (KLU)                        902                                   0.233                         0.107
NR (KLU)                               885                                   0.25                          0.123
FDPF XB (SLU)                          808                                   0.384                         0.267
FDPF BX (SLU)                          770                                   0.419                         0.297
FDPF XB (KLU)                          838                                   0.324                         0.205
FDPF BX (KLU)                          805                                   0.348                         0.226
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
FDPF XB (SLU)                             6.1e-05        1.91e-06          1.53e-05
FDPF BX (SLU)                             6.1e-05        1.91e-06          7.63e-06
FDPF XB (KLU)                             6.1e-05        1.91e-06          1.53e-05
FDPF BX (KLU)                             6.1e-05        1.91e-06          7.63e-06
=================================  ==============  ==============  ================