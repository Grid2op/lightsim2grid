Lightsim2grid 0.8.1 and grid2op 1.10.0 (python 3.11)
====================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-03-26 11:42  CET
- system: Linux 5.4.0-153-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-5960X CPU @ 3.00GHz
- python version: 3.11.8.final.0 (64 bit)
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
PP                                      32.9                               27.4                         15.8
GS                                     729                                  0.499                        0.396
GS synch                               728                                  0.507                        0.404
NR single (SLU)                        912                                  0.196                        0.0844
NR (SLU)                               935                                  0.193                        0.0846
NR single (KLU)                       1000                                  0.128                        0.0232
NR (KLU)                               989                                  0.128                        0.0218
NR single (NICSLU *)                   996                                  0.128                        0.0232
NR (NICSLU *)                         1000                                  0.126                        0.0212
NR single (CKTSO *)                   1010                                  0.126                        0.0217
NR (CKTSO *)                          1010                                  0.124                        0.0199
FDPF XB (SLU)                          992                                  0.138                        0.0347
FDPF BX (SLU)                          974                                  0.153                        0.0492
FDPF XB (KLU)                          996                                  0.132                        0.0286
FDPF BX (KLU)                          988                                  0.144                        0.0414
FDPF XB (NICSLU *)                     978                                  0.135                        0.0301
FDPF BX (NICSLU *)                     973                                  0.146                        0.0426
FDPF XB (CKTSO *)                      994                                  0.132                        0.0291
FDPF BX (CKTSO *)                      989                                  0.144                        0.0412
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

- date: 2024-03-26 11:49  CET
- system: Linux 5.4.0-153-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-5960X CPU @ 3.00GHz
- python version: 3.11.8.final.0 (64 bit)
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
PP                                      29.1                                31.2                          19.3
GS                                       4.12                              241                           241
GS synch                                29.6                                32.7                          32.6
NR single (SLU)                        531                                   0.888                         0.749
NR (SLU)                               513                                   0.929                         0.784
NR single (KLU)                        785                                   0.284                         0.148
NR (KLU)                               793                                   0.268                         0.133
NR single (NICSLU *)                   784                                   0.274                         0.137
NR (NICSLU *)                          795                                   0.262                         0.125
NR single (CKTSO *)                    784                                   0.27                          0.132
NR (CKTSO *)                           812                                   0.25                          0.116
FDPF XB (SLU)                          728                                   0.386                         0.256
FDPF BX (SLU)                          714                                   0.409                         0.279
FDPF XB (KLU)                          750                                   0.345                         0.215
FDPF BX (KLU)                          737                                   0.367                         0.236
FDPF XB (NICSLU *)                     724                                   0.354                         0.219
FDPF BX (NICSLU *)                     740                                   0.361                         0.231
FDPF XB (CKTSO *)                      756                                   0.339                         0.209
FDPF BX (CKTSO *)                      744                                   0.359                         0.229
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
