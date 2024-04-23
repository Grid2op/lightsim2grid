Lightsim2grid 0.8.2 and grid2op 1.10.1 (python 3.8)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 11:10  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.10.1
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
PP                                      36.6                               24.6                         17.2
GS                                     819                                  0.432                        0.339
GS synch                               809                                  0.445                        0.353
NR single (SLU)                       1040                                  0.167                        0.0686
NR (SLU)                              1040                                  0.167                        0.0694
NR single (KLU)                       1110                                  0.114                        0.0197
NR (KLU)                              1100                                  0.113                        0.0185
NR single (NICSLU *)                  1120                                  0.112                        0.0194
NR (NICSLU *)                         1130                                  0.11                         0.0175
NR single (CKTSO *)                   1180                                  0.105                        0.0171
NR (CKTSO *)                          1130                                  0.109                        0.0165
FDPF XB (SLU)                         1080                                  0.124                        0.0294
FDPF BX (SLU)                         1080                                  0.135                        0.0419
FDPF XB (KLU)                         1100                                  0.116                        0.024
FDPF BX (KLU)                         1090                                  0.127                        0.0345
FDPF XB (NICSLU *)                    1120                                  0.117                        0.0245
FDPF BX (NICSLU *)                    1090                                  0.128                        0.0355
FDPF XB (CKTSO *)                     1180                                  0.108                        0.0224
FDPF BX (CKTSO *)                     1160                                  0.119                        0.0322
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

- date: 2024-04-23 11:16  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
- numpy version: 1.24.4
- pandas version: 2.0.3
- pandapower version: 2.14.6
- grid2op version: 1.10.1
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
PP                                      32.9                                27.5                         19.7
GS                                       4.63                              215                          215
GS synch                                34.8                                27.9                         27.7
NR single (SLU)                        582                                   0.793                        0.666
NR (SLU)                               600                                   0.778                        0.654
NR single (KLU)                        950                                   0.227                        0.116
NR (KLU)                               963                                   0.215                        0.103
NR single (NICSLU *)                   961                                   0.22                         0.109
NR (NICSLU *)                          969                                   0.209                        0.0976
NR single (CKTSO *)                    967                                   0.214                        0.104
NR (CKTSO *)                           975                                   0.204                        0.0919
FDPF XB (SLU)                          877                                   0.315                        0.207
FDPF BX (SLU)                          866                                   0.333                        0.225
FDPF XB (KLU)                          912                                   0.276                        0.168
FDPF BX (KLU)                          905                                   0.29                         0.184
FDPF XB (NICSLU *)                     908                                   0.279                        0.171
FDPF BX (NICSLU *)                     899                                   0.291                        0.185
FDPF XB (CKTSO *)                      910                                   0.275                        0.168
FDPF BX (CKTSO *)                      893                                   0.293                        0.185
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
