Lightsim2grid 0.8.1 and grid2op 1.10.1
======================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-03-26 10:52  CET
- system: Linux 5.4.0-153-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-5960X CPU @ 3.00GHz
- python version: 3.10.13.final.0 (64 bit)
- numpy version: 1.23.5
- pandas version: 2.2.1
- pandapower version: 2.13.1
- grid2op version: 1.10.1
- lightsim2grid version: 0.8.1
- lightsim2grid extra information: 

	- klu_solver_available: True 
	- nicslu_solver_available: True 
	- cktso_solver_available: True 
	- compiled_march_native: True 
	- compiled_o3_optim: True 

====================  ======================  ===================================  ============================
case14_sandbox          grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
====================  ======================  ===================================  ============================
PP                                      41.9                               20.5                          7.14
GS                                     678                                  0.507                        0.395
GS synch                               677                                  0.506                        0.393
NR single (SLU)                        858                                  0.198                        0.0823
NR (SLU)                               832                                  0.208                        0.0856
NR single (KLU)                        910                                  0.136                        0.0232
NR (KLU)                               909                                  0.134                        0.0213
NR single (NICSLU *)                   898                                  0.138                        0.023
NR (NICSLU *)                          898                                  0.135                        0.0212
NR single (CKTSO *)                    895                                  0.137                        0.0221
NR (CKTSO *)                           903                                  0.134                        0.0198
FDPF XB (SLU)                          894                                  0.148                        0.0351
FDPF BX (SLU)                          890                                  0.16                         0.0489
FDPF XB (KLU)                          885                                  0.143                        0.0294
FDPF BX (KLU)                          891                                  0.151                        0.0412
FDPF XB (NICSLU *)                     891                                  0.143                        0.0299
FDPF BX (NICSLU *)                     881                                  0.156                        0.0425
FDPF XB (CKTSO *)                      890                                  0.142                        0.0295
FDPF BX (CKTSO *)                      892                                  0.152                        0.0414
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

- date: 2024-03-26 10:59  CET
- system: Linux 5.4.0-153-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-5960X CPU @ 3.00GHz
- python version: 3.10.13.final.0 (64 bit)
- numpy version: 1.23.5
- pandas version: 2.2.1
- pandapower version: 2.13.1
- grid2op version: 1.10.1
- lightsim2grid version: 0.8.1
- lightsim2grid extra information: 

	- klu_solver_available: True 
	- nicslu_solver_available: True 
	- cktso_solver_available: True 
	- compiled_march_native: True 
	- compiled_o3_optim: True 

=====================  ======================  ===================================  ============================
neurips_2020_track2      grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
=====================  ======================  ===================================  ============================
PP                                      38.3                                22.5                           9.18
GS                                       3.81                              261                           261
GS synch                                29.2                                33.2                          33
NR single (SLU)                        494                                   0.912                         0.76
NR (SLU)                               484                                   0.935                         0.777
NR single (KLU)                        709                                   0.302                         0.155
NR (KLU)                               714                                   0.29                          0.143
NR single (NICSLU *)                   706                                   0.295                         0.146
NR (NICSLU *)                          717                                   0.28                          0.132
NR single (CKTSO *)                    712                                   0.286                         0.138
NR (CKTSO *)                           722                                   0.272                         0.126
FDPF XB (SLU)                          658                                   0.414                         0.271
FDPF BX (SLU)                          645                                   0.436                         0.292
FDPF XB (KLU)                          680                                   0.366                         0.225
FDPF BX (KLU)                          667                                   0.388                         0.244
FDPF XB (NICSLU *)                     674                                   0.364                         0.22
FDPF BX (NICSLU *)                     665                                   0.384                         0.239
FDPF XB (CKTSO *)                      677                                   0.362                         0.22
FDPF BX (CKTSO *)                      665                                   0.383                         0.24
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

