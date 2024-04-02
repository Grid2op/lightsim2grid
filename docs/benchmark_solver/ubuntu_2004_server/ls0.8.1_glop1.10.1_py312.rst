Lightsim2grid 0.8.1 and grid2op 1.10.1 (python 3.12)
====================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-03-26 13:42  CET
- system: Linux 5.4.0-153-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-5960X CPU @ 3.00GHz
- python version: 3.12.2.final.0 (64 bit)
- numpy version: 1.26.4
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
PP                                      31.3                               28.9                         16.8
GS                                     716                                  0.502                        0.393
GS synch                               696                                  0.522                        0.41
NR single (SLU)                        912                                  0.196                        0.0829
NR (SLU)                               911                                  0.197                        0.084
NR single (KLU)                        971                                  0.133                        0.023
NR (KLU)                               973                                  0.132                        0.0213
NR single (NICSLU *)                   955                                  0.136                        0.023
NR (NICSLU *)                          952                                  0.135                        0.0213
NR single (CKTSO *)                    970                                  0.132                        0.0216
NR (CKTSO *)                           975                                  0.13                         0.0199
FDPF XB (SLU)                          961                                  0.145                        0.0355
FDPF BX (SLU)                          949                                  0.158                        0.0491
FDPF XB (KLU)                          968                                  0.137                        0.0285
FDPF BX (KLU)                          957                                  0.15                         0.0412
FDPF XB (NICSLU *)                     963                                  0.139                        0.0295
FDPF BX (NICSLU *)                     943                                  0.154                        0.0429
FDPF XB (CKTSO *)                      971                                  0.138                        0.0291
FDPF BX (CKTSO *)                      956                                  0.151                        0.0416
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

- date: 2024-03-26 13:51  CET
- system: Linux 5.4.0-153-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-5960X CPU @ 3.00GHz
- python version: 3.12.2.final.0 (64 bit)
- numpy version: 1.26.4
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
PP                                      28.1                                32.3                          20.4
GS                                       4.05                              245                           245
GS synch                                29.4                                33                            32.9
NR single (SLU)                        519                                   0.9                           0.754
NR (SLU)                               518                                   0.917                         0.771
NR single (KLU)                        763                                   0.293                         0.151
NR (KLU)                               770                                   0.28                          0.139
NR single (NICSLU *)                   777                                   0.279                         0.139
NR (NICSLU *)                          774                                   0.27                          0.127
NR single (CKTSO *)                    792                                   0.27                          0.13
NR (CKTSO *)                           786                                   0.26                          0.119
FDPF XB (SLU)                          700                                   0.406                         0.268
FDPF BX (SLU)                          688                                   0.431                         0.294
FDPF XB (KLU)                          725                                   0.36                          0.223
FDPF BX (KLU)                          719                                   0.376                         0.24
FDPF XB (NICSLU *)                     730                                   0.351                         0.215
FDPF BX (NICSLU *)                     732                                   0.365                         0.229
FDPF XB (CKTSO *)                      748                                   0.343                         0.209
FDPF BX (CKTSO *)                      729                                   0.367                         0.232
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

