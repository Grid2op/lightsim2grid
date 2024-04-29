Lightsim2grid 0.8.2 and grid2op 1.9.7 (python 3.9)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 12:06  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.9.19.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.7
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
PP                                      34.9                               22.1                         14.4
GS                                     836                                  0.433                        0.337
GS synch                               827                                  0.44                         0.346
NR single (SLU)                       1070                                  0.167                        0.0691
NR (SLU)                              1070                                  0.168                        0.0702
NR single (KLU)                       1140                                  0.114                        0.0195
NR (KLU)                              1130                                  0.113                        0.0181
NR single (NICSLU *)                  1150                                  0.112                        0.0187
NR (NICSLU *)                         1160                                  0.11                         0.0172
NR single (CKTSO *)                   1090                                  0.118                        0.019
NR (CKTSO *)                          1090                                  0.116                        0.0175
FDPF XB (SLU)                         1130                                  0.122                        0.0293
FDPF BX (SLU)                         1120                                  0.135                        0.0419
FDPF XB (KLU)                         1140                                  0.117                        0.0243
FDPF BX (KLU)                         1120                                  0.129                        0.0352
FDPF XB (NICSLU *)                    1140                                  0.118                        0.0247
FDPF BX (NICSLU *)                    1100                                  0.131                        0.0364
FDPF XB (CKTSO *)                     1080                                  0.124                        0.0257
FDPF BX (CKTSO *)                     1070                                  0.134                        0.0366
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

- date: 2024-04-23 12:12  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.9.19.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.9.7
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
PP                                      15.7                                26.5                         17.4
GS                                       4.59                              217                          217
GS synch                                34.9                                27.7                         27.6
NR single (SLU)                        616                                   0.756                        0.633
NR (SLU)                               580                                   0.86                         0.65
NR single (KLU)                        904                                   0.243                        0.124
NR (KLU)                               915                                   0.23                         0.111
NR single (NICSLU *)                   910                                   0.234                        0.115
NR (NICSLU *)                          934                                   0.221                        0.103
NR single (CKTSO *)                    977                                   0.214                        0.102
NR (CKTSO *)                           978                                   0.206                        0.0928
FDPF XB (SLU)                          833                                   0.337                        0.222
FDPF BX (SLU)                          822                                   0.357                        0.242
FDPF XB (KLU)                          863                                   0.298                        0.183
FDPF BX (KLU)                          850                                   0.315                        0.2
FDPF XB (NICSLU *)                     910                                   0.283                        0.173
FDPF BX (NICSLU *)                     904                                   0.296                        0.188
FDPF XB (CKTSO *)                      914                                   0.281                        0.172
FDPF BX (CKTSO *)                      902                                   0.296                        0.187
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
