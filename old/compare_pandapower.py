# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of PyKLU2Grid, PyKLU2Grid a implements a c++ backend targeting the Grid2Op platform.

import numpy as np
import sys
from scipy import sparse
from lightsim2grid.compute_powerflow import KLUSolver

import pdb
for it in range(1,5):
    Ybus = np.load("Ybus_{}.npy".format(it))
    Ybus = sparse.csc_matrix(Ybus)
    V = np.load("V_{}.npy".format(it))
    dS_dVm_pp_ = np.load("dS_dVm_{}.npy".format(it))
    dS_dVm_pp = sparse.csc_matrix(dS_dVm_pp_)
    dS_dVa_pp_ = np.load("dS_dVa_{}.npy".format(it))
    dS_dVa_pp = sparse.csc_matrix(dS_dVa_pp_)

    solver = KLUSolver()
    dS_dVm, dS_dVa = solver.get_ds_test(Ybus, V)

    comp_val_Vm = np.abs(dS_dVm - dS_dVm_pp)
    comp_val_Va = np.abs(dS_dVa - dS_dVa_pp)
    print("Results for iteration {}".format(it))
    print("\tMaximum difference (absolute value) for dS_dVm: {}".format(np.max(np.abs(comp_val_Vm))))
    print("\tMaximum difference (absolute value) for dS_dVa: {}".format(np.max(np.abs(comp_val_Va))))
pdb.set_trace()

sys.exit()
# test the full newton raphson (use test_pyklu instead now)
it_num = 0

solver = KLUSolver()
solver2 = KLUSolver()
solver3 = KLUSolver()
# test the Newton Raphson solver
tol = 1e-08
Ybus = np.load("Ybus.npy")
Ybus = sparse.csc_matrix(Ybus)
Sbus = np.load("Sbus.npy")
V0 = np.load("V0.npy")
pv = np.load("pv.npy")
pq = np.load("pq.npy")
pvpq = np.r_[pv, pq]
max_it = 10
tol = 1e-08

V = 1.0 * V0
F = np.zeros(181)
has_conv = solver.do_newton(Ybus, V, Sbus, pv, pq, max_it, tol)
error_status = solver.get_error()
J = solver.get_J()
J_pp_ = np.load("J_{}.npy".format(4))
# F_pp = np.load("F_{}.npy".format(4))
J_pp = sparse.csc_matrix(J_pp_)
# has_conv_pandapower = solver2.initialize_test(J_pp)
print("Has the solver converged: {}".format(has_conv))
print("In how many iteration: {}".format(solver.get_nb_iter()))

test2 = np.where(J_pp.toarray() != 0)
test = np.where(J.toarray() != 0)

print("Are the non null values identical to what they should be?")
print("\t for rows: {}".format(np.all(test[0] == test2[0])))
print("\t for columns: {}".format(np.all(test[1] == test2[1])))

comp_val = np.abs(J - J_pp)
comp_val = comp_val.toarray()
print("Is the jacobian the same?")
print("\t for J11 (dS_dVa_r): {}".format(np.sum(np.abs(comp_val[:len(pvpq), :len(pvpq)]))))
print("\t for J21 (dS_dVa_i): {}".format(np.sum(np.abs(comp_val[len(pvpq):, :len(pvpq)]))))
print("\t for J12 (dS_dVm_r): {}".format(np.sum(np.abs(comp_val[:len(pvpq), len(pvpq):]))))
print("\t for J22 (dS_dVm_i): {}".format(np.sum(np.abs(comp_val[len(pvpq):, len(pvpq):]))))
sys.exit()

is_init = False
# check the "create jacobian" stuff
tol = 1e-08
Ybus = np.load("Ybus.npy")
Ybus = sparse.csc_matrix(Ybus)
Sbus = np.load("Sbus.npy")
V0 = np.load("V0.npy")
pv = np.load("pv.npy")
pq = np.load("pq.npy")

# other variable initialized from above
iwamoto = False
numba = False
max_it = 10
i = 0
V = V0
Va = np.angle(V)
Vm = abs(V)
pvpq = np.r_[pv, pq]
pvpq_lookup = np.zeros(max(Ybus.indices) + 1, dtype=int)
pvpq_lookup[pvpq] = np.arange(len(pvpq))
npv = len(pv)
npq = len(pq)
j1 = 0
j2 = npv  # j1:j2 - V angle of pv buses
j3 = j2
j4 = j2 + npq  # j3:j4 - V angle of pq buses
j5 = j4
j6 = j4 + npq  # j5:j6 - V mag of pq buses
createJ = get_fastest_jacobian_function(pvpq, pq, numba)

# mimic the code
F = solver._evaluate_Fx(Ybus, V, Sbus, pv, pq)
converged = solver._check_for_convergence(F, tol)
while (not converged and i < max_it):
    i = i + 1
    J = solver.create_jacobian_matrix(Ybus, V, pq, pvpq)
    # to test
    J2_ = create_jacobian_matrix(Ybus, V, pvpq, pq, createJ, pvpq_lookup, npv, npq, numba)
    J2 = sparse.csc_matrix(J2_)

    J_pp_ = np.load("J_{}.npy".format(i))
    J_pp = sparse.csc_matrix(J_pp_)

    test2 = np.where(J2.toarray() != 0)
    test = np.where(J.toarray() != 0)

    print("Are the non null values identical: ")
    print("\t for rows: {}".format(np.all(test[0] == test2[0])))
    print("\t for columns: {}".format(np.all(test[1] == test2[1])))


    comp_val = np.abs(J - J2)
    comp_val = comp_val.toarray()
    print("Is J the same:")
    print("\t for J11 (dS_dVa_r): {}".format(np.sum(np.abs(comp_val[:len(pvpq), :len(pvpq)]))))
    print("\t for J21 (dS_dVa_i): {}".format(np.sum(np.abs(comp_val[len(pvpq):, :len(pvpq)]))))
    print("\t for J12 (dS_dVm_r): {}".format(np.sum(np.abs(comp_val[:len(pvpq), len(pvpq):]))))
    print("\t for J22 (dS_dVm_i): {}".format(np.sum(np.abs(comp_val[len(pvpq):, len(pvpq):]))))
    pdb.set_trace()
    # J = sparse.csc_matrix(J)

sys.exit()

# check the one_iter function
tol = 1e-08
Ybus = np.load("Ybus.npy")
Ybus = sparse.csc_matrix(Ybus)
Sbus = np.load("Sbus.npy")
V0 = np.load("V0.npy")
pv = np.load("pv.npy")
pq = np.load("pq.npy")

# other variable initialized from above
iwamoto = False
numba = False
max_it = 10
i = 0
V = V0
Va = np.angle(V)
Vm = abs(V)
pvpq = np.r_[pv, pq]
pvpq_lookup = np.zeros(max(Ybus.indices) + 1, dtype=int)
pvpq_lookup[pvpq] = np.arange(len(pvpq))
npv = len(pv)
npq = len(pq)
j1 = 0
j2 = npv  # j1:j2 - V angle of pv buses
j3 = j2
j4 = j2 + npq  # j3:j4 - V angle of pq buses
j5 = j4
j6 = j4 + npq  # j5:j6 - V mag of pq buses
createJ = get_fastest_jacobian_function(pvpq, pq, numba)

# mimic the code
F = solver._evaluate_Fx(Ybus, V, Sbus, pv, pq)
converged = solver._check_for_convergence(F, tol)
while (not converged and i < max_it):
    i = i + 1
    J = create_jacobian_matrix(Ybus, V, pvpq, pq, createJ, pvpq_lookup, npv, npq, numba)
    J = sparse.csc_matrix(J)
    if not is_init:
        is_init = True
        solver.analyze(J)

    J_pp = np.load("J_{}.npy".format(i))
    F_pp = np.load("F_{}.npy".format(i))
    print("assertion is equal to pandapower for iteration {}".format(i))
    print("\tJ: {}".format(np.all(np.abs(J - J_pp) <= 1e-5)))
    same_as_pp = np.all(np.abs(F - F_pp) <= 1e-5)
    print("\tF: {}".format(same_as_pp))
    if not same_as_pp:
        pdb.set_trace()

    F, V = solver.one_iter(J, F, pv, pq,
                           V, #Va.astype(np.float), Vm.astype(np.float),
                           Ybus, Sbus)
    converged = solver._check_for_convergence(F, tol)

sys.exit()


# check evaluate function and norm
tol = 1e-08
Ybus = np.load("Ybus.npy")
Ybus = sparse.csc_matrix(Ybus)
Sbus = np.load("Sbus.npy")
V0 = np.load("V0.npy")
pv = np.load("pv.npy")
pq = np.load("pq.npy")

# other variable initialized from above
iwamoto = False
numba = False
max_it = 10
i = 0
V = V0
Va = np.angle(V)
Vm = abs(V)
pvpq = np.r_[pv, pq]
pvpq_lookup = np.zeros(max(Ybus.indices) + 1, dtype=int)
pvpq_lookup[pvpq] = np.arange(len(pvpq))
npv = len(pv)
npq = len(pq)
j1 = 0
j2 = npv  # j1:j2 - V angle of pv buses
j3 = j2
j4 = j2 + npq  # j3:j4 - V angle of pq buses
j5 = j4
j6 = j4 + npq  # j5:j6 - V mag of pq buses
createJ = get_fastest_jacobian_function(pvpq, pq, numba)

# mimic the code
F = solver._evaluate_Fx(Ybus, V, Sbus, pv, pq)
converged = solver._check_for_convergence(F, tol)
while (not converged and i < max_it):
    i = i + 1
    J = create_jacobian_matrix(Ybus, V, pvpq, pq, createJ, pvpq_lookup, npv, npq, numba)
    J = sparse.csc_matrix(J)
    if not is_init:
        is_init = True
        solver.analyze(J)

    # solve the system
    dx = 1.0 * F
    solver.solve(J, dx)
    dx *= -1.0

    # update voltage
    if npv and not iwamoto:
        Va[pv] = Va[pv] + dx[j1:j2]
    if npq and not iwamoto:
        Va[pq] = Va[pq] + dx[j3:j4]
        Vm[pq] = Vm[pq] + dx[j5:j6]

    V = Vm * np.exp(1j * Va)
    Vm = abs(V)  # update Vm and Va again in case
    Va = np.angle(V)  # we wrapped around with a negative Vm

    J_pp = np.load("J_{}.npy".format(i))
    dx_pp = np.load("dx_{}.npy".format(i))
    F_pp = np.load("F_{}.npy".format(i))

    print("assertion is equal to pandapower for iteration {}".format(i))
    print("\tJ: {}".format(np.all(np.abs(J - J_pp) <= 1e-5)))
    print("\tdx: {}".format(np.all(np.abs(dx - dx_pp) <= 1e-5)))
    print("\tF: {}".format(np.all(np.abs(F - F_pp) <= 1e-5)))

    F = solver._evaluate_Fx(Ybus, V, Sbus, pv, pq)
    converged = solver._check_for_convergence(F, tol)




sys.exit()
# check the invert solver
for it_num in [1, 2, 3, 4]:
    J = np.load("J_{}.npy".format(it_num))
    dx = np.load("dx_{}.npy".format(it_num))
    F = np.load("F_{}.npy".format(it_num))

    F_klu = 1.0 * F

    # need to be in csc matrix
    compress_A = sparse.csc_matrix(J)
    if not is_init:
        is_init = True
        # Ap = compress_A.indptr
        # Ai = compress_A.indices
        # n = Ap.size - 1
        # solver.analyze(int(n), Ap, Ai)
        solver.analyze(compress_A)  #int(n), Ap, Ai)

    # print(F_klu[:5])
    ## solver.solve(Ap, Ai, compress_A.data, F_klu)
    solver.solve(compress_A, F_klu)
    dx_klu = -1.0 * F_klu
    # print(F_klu[:5])

    # res_klu = np.matmul(J, -1.0 * dx_klu) - F
    # res = np.matmul(J, -1.0 * dx) - F
    # print("res_klu: {}".format(np.sum(np.abs(res_klu))))
    # print("res: {}".format(np.sum(np.abs(res))))
    # print(F_klu[:5])
sys.exit()

# # dx = -1 * spsolve(J, F, permc_spec=permc_spec, use_umfpack=use_umfpack)
# see pandapower.pypower.newtonpf and newtonpf
# import numpy as np
# F = np.array([-2.40868988e+01,  1.28470484e+01, -5.90647951e+00, -2.84961091e+01,
#                 3.43054654e-01, -1.50490473e+01, -2.47184985e+00, -1.79735240e+00,
#                 1.33656911e+01,  7.72235439e+00,  3.37569963e+01, -4.73176627e+00,
#                 2.70588869e+00, -1.19207512e+01,  3.79737055e+00,  1.98062446e+01,
#                 8.28495595e+00, -1.81921518e+01,  2.51794141e+01, -1.72386939e+01,
#                -1.37021142e+01, -4.74929487e+00, -2.00634711e+01, -1.54528301e+01,
#                 1.14227448e+00,  1.00922067e+01,  4.37165896e+01,  4.16234337e+00,
#                -1.04932999e+01,  1.04290700e+00, -1.50213339e+01,  7.27652222e+00,
#                -2.59131600e+00,  3.88186912e+00,  8.74638448e+00,  4.15407621e+01,
#                -6.71110016e+00,  1.66592934e+01, -1.61290061e+00,  8.06365233e-01,
#                -1.68274547e+01, -2.37032325e+01,  1.73276739e+01,  6.94067069e+01,
#                -1.47350705e+01,  1.59962144e+00, -4.44352738e+00,  8.02846245e+00,
#                 3.13932398e+01, -2.63415646e+00, -7.43639078e+00, -4.24722772e+00,
#                 2.69336655e+01,  1.12731147e+01,  4.08917007e+00, -2.36993936e+00,
#                 5.27562435e-01,  1.58100217e+01,  1.23632415e+01,  7.94204158e+00,
#                 9.86388471e+00,  8.07155493e+00,  1.45914180e+00,  2.68449448e+01,
#                 5.84277115e+00,  4.10618936e+00,  1.98248611e+00,  3.85203561e+00,
#                 5.79307800e+00,  2.35432909e+01,  1.15023935e+00,  1.30073890e+00,
#                -2.05860715e+01,  7.13510424e+00,  2.24722122e+01, -9.18612791e+00,
#                 6.96005994e+00,  2.21490266e+01,  7.16040052e-01, -1.80440712e-01,
#                 1.17984140e+01, -4.52227074e+01,  1.62415142e+01,  2.25400037e+00,
#                 1.38603338e+00, -2.21561568e+00, -2.05808708e+01, -2.18683678e+01,
#                -1.34011415e+01, -6.56817916e+00,  2.81057920e+00,  7.88158904e+00,
#                 1.27255644e+01,  1.29745451e+01, -5.83253646e+00,  7.17083252e-01,
#                 5.09162251e-01, -5.94932087e+00, -1.43211727e+01, -4.04189491e+01,
#                 4.72553476e+00,  1.23653555e+01, -3.02116063e+01, -1.51582805e+01,
#                -3.13453808e+01, -3.21527641e+00, -1.64994144e-01,  6.66301252e+00,
#                 5.61436619e+00,  1.17732094e+00, -1.10505433e+00,  1.84481045e+00,
#                -8.76751187e+00,  4.44977293e+00, -2.71642610e+00, -1.01924691e+01,
#                -1.18987380e+01,  4.41585924e+01,  3.09139394e+01, -3.01378151e+00,
#                 2.84169429e+00,  7.44351034e+01,  3.11849237e+01,  2.36163415e+01,
#                 4.50154050e+01,  3.82588132e+01,  1.07337802e+01,  1.01060329e+02,
#                 2.90804711e+01,  1.05846354e+01,  1.25022166e+01, -7.41405972e+01,
#                 2.70390717e+01,  8.95997795e+01,  9.71550043e+00,  6.70690990e+00,
#                -5.26990626e+01,  3.23719114e+01,  7.01225962e+01, -2.73649716e+00,
#                 2.63237860e+01,  1.06000375e+02, -1.44978425e+02,  1.53015617e+02,
#                 4.77992116e+01, -2.65001576e+02,  5.80383447e+01,  9.50382786e+00,
#                -1.61205350e+00,  1.99713517e+00, -6.89789192e+01, -6.78410937e+01,
#                -3.23006253e+01, -8.97726650e+00,  8.78268413e+00,  3.94574468e+01,
#                 3.75451783e+01,  3.73943670e+01, -3.17081714e+01,  3.45927260e+00,
#                 1.31512848e+02, -3.08275863e+01, -5.04497872e+01, -4.77448839e+02,
#                 1.49139959e+01,  5.79062175e+01, -7.54945566e+01, -4.57182853e+01,
#                -3.16657194e+02,  9.59132472e+00, -3.02318224e+00,  9.00065162e+00,
#                 1.62118985e+01,  3.27317969e+00,  6.00983438e+00,  1.02787563e+01,
#                -2.23377791e+01,  4.15849366e+01, -1.01716755e+01, -4.14302504e+01,
#                -4.62477040e+01])
#
# dx = np.array([-0.05463578, -0.0382218 , -0.0323992 , -0.03340859, -0.03854723,
#                -0.04248839, -0.05978398, -0.040956  , -0.05303533, -0.05499707,
#                -0.01695892, -0.05305905, -0.0530527 , -0.05076724, -0.03370322,
#                -0.05496989, -0.05634815, -0.04660005, -0.0575151 , -0.05110835,
#                -0.04559774, -0.04825972, -0.04809417, -0.04924428, -0.04541869,
#                -0.02684788, -0.02940566, -0.05857076, -0.02496616, -0.02521909,
#                -0.02509464, -0.02507503, -0.02213192, -0.02275906, -0.02353401,
#                -0.03105105, -0.00906439, -0.05941174, -0.02118386, -0.0130712 ,
#                -0.00438955, -0.00364688, -0.01348649, -0.02349328, -0.01619588,
#                -0.02293209, -0.01981847, -0.02372939, -0.0842955 , -0.02124834,
#                -0.0212674 , -0.02777709, -0.03304965, -0.05715996, -0.05782893,
#                -0.0268229 , -0.02310438, -0.03333692, -0.03722006, -0.03864532,
#                -0.04614045, -0.04621666, -0.05869014, -0.00326461, -0.05548527,
#                -0.05853158, -0.05774513, -0.05542647, -0.04826335, -0.05596753,
#                -0.04529825, -0.04251601, -0.04086443, -0.04841272, -0.05051452,
#                -0.05307107, -0.0508482 , -0.04820202, -0.04918101, -0.04491326,
#                -0.04831426, -0.05807827, -0.04809035, -0.04276021, -0.03433976,
#                -0.02963108, -0.02460523, -0.02941373, -0.0266234 , -0.02437822,
#                -0.02450553, -0.02503458, -0.02529449, -0.02485161, -0.02247777,
#                -0.05901855, -0.02266533, -0.02235809, -0.02510739, -0.01683595,
#                -0.01211087, -0.00372743, -0.0134815 , -0.01519699, -0.07082179,
#                -0.01870511, -0.0190566 , -0.01928269, -0.01627639, -0.01918528,
#                -0.01742291, -0.02224974, -0.02401617, -0.02099645, -0.02165357,
#                -0.02148729, -0.02515766, -0.01806002, -0.00459003,  0.00165252,
#                -0.0006491 , -0.02782071, -0.02303107, -0.02223859, -0.02883201,
#                -0.02889966, -0.01571272, -0.03883342, -0.02101821, -0.00616555,
#                -0.00580896,  0.00540714, -0.03163609, -0.02150471, -0.03076257,
#                -0.01973142,  0.0099948 , -0.02745302, -0.02587349, -0.00412821,
#                -0.01837788, -0.0092111 ,  0.00112257, -0.02734521, -0.01934987,
#                 0.0123894 , -0.02243306, -0.01239446, -0.0051167 , -0.00324492,
#                 0.02818471,  0.03183937,  0.01174089, -0.02200706, -0.03166426,
#                -0.04196828, -0.01856925, -0.02970183,  0.00337246, -0.00046795,
#                -0.02056048, -0.00585641,  0.03083093,  0.01365736, -0.00293157,
#                -0.02171797,  0.01385313,  0.01979025,  0.05623961,  0.00717467,
#                -0.00108937, -0.00528149, -0.00996533, -0.00306131, -0.00233939,
#                -0.00424623,  0.00021747, -0.00913819,  0.0027597 ,  0.02196277,
#                 0.03493246])

import lightsim2grid
dx_klu = -1.0 * lightsim2grid.solve_linear_system(J, F_klu)
# from scipy import sparse
# import ctypes
# libklu = ctypes.cdll.LoadLibrary('libpyklu.so')
#
# compress_A = sparse.csc_matrix(J)
# Ap = compress_A.indptr
# Ai = compress_A.indices
# Ax = compress_A.data
# n = Ap.size - 1
# c_Ap = np.ctypeslib.as_ctypes(Ap)
# c_Ai = np.ctypeslib.as_ctypes(Ai)
# c_Ax = np.ctypeslib.as_ctypes(Ax)
# c_b = np.ctypeslib.as_ctypes(F_klu)
# libklu.solve_linear_system(n, c_Ap, c_Ai, c_Ax, c_b)
# dx_klu = -1.0 * np.array(c_b)


res_klu = np.matmul(J, -1.0 * dx_klu) - F
res_pp = np.matmul(J, -1.0 * dx) - F

import pdb
pdb.set_trace()
# with open("F.npy", "r") as f:
#     F = np.load(f.read())
# with open("J.npy", "r") as f:
#     J = np.load(f)
# with open("dx.npy", "r") as f:
#     dx = np.load(f)