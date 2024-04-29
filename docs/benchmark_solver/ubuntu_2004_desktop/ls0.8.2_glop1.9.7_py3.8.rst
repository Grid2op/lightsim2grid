Lightsim2grid 0.8.2 and grid2op 1.9.7 (python 3.8)
=================================================================

l2rpn_case14_sandbox
---------------------

Configuration:

- date: 2024-04-23 10:50  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
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
PP                                      32.3                               24.4                         16.7
GS                                     822                                  0.426                        0.331
GS synch                               794                                  0.446                        0.349
NR single (SLU)                       1010                                  0.169                        0.0683
NR (SLU)                              1050                                  0.164                        0.0665
NR single (KLU)                       1080                                  0.117                        0.0195
NR (KLU)                              1070                                  0.117                        0.018
NR single (NICSLU *)                  1010                                  0.125                        0.0203
NR (NICSLU *)                         1040                                  0.119                        0.0179
NR single (CKTSO *)                   1130                                  0.111                        0.0174
NR (CKTSO *)                          1140                                  0.109                        0.0159
FDPF XB (SLU)                         1070                                  0.126                        0.029
FDPF BX (SLU)                         1050                                  0.14                         0.041
FDPF XB (KLU)                         1080                                  0.12                         0.0239
FDPF BX (KLU)                         1060                                  0.131                        0.034
FDPF XB (NICSLU *)                    1130                                  0.115                        0.0229
FDPF BX (NICSLU *)                    1110                                  0.127                        0.033
FDPF XB (CKTSO *)                     1130                                  0.114                        0.0226
FDPF BX (CKTSO *)                     1120                                  0.124                        0.0321
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

- date: 2024-04-23 10:56  CEST
- system: Linux 5.15.0-105-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
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
PP                                      14.9                                29.5                          20.2
GS                                       4.26                              234                           233
GS synch                                34.6                                28                            27.9
NR single (SLU)                        605                                   0.759                         0.633
NR (SLU)                               596                                   0.78                          0.652
NR single (KLU)                        878                                   0.245                         0.124
NR (KLU)                               891                                   0.233                         0.112
NR single (NICSLU *)                   937                                   0.226                         0.111
NR (NICSLU *)                          949                                   0.211                         0.098
NR single (CKTSO *)                    942                                   0.217                         0.103
NR (CKTSO *)                           954                                   0.206                         0.091
FDPF XB (SLU)                          819                                   0.331                         0.215
FDPF BX (SLU)                          803                                   0.35                          0.233
FDPF XB (KLU)                          845                                   0.291                         0.175
FDPF BX (KLU)                          829                                   0.31                          0.193
FDPF XB (NICSLU *)                     890                                   0.275                         0.164
FDPF BX (NICSLU *)                     884                                   0.29                          0.18
FDPF XB (CKTSO *)                      902                                   0.274                         0.165
FDPF BX (CKTSO *)                      883                                   0.289                         0.179
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
