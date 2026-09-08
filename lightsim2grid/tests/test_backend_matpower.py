# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
`LightSimBackend(loader_method="matpower")`: a grid2op environment whose powergrid is
a MATPOWER case (a `grid.m` / `grid.mat` file), rather than a pandapower `grid.json`
or a pypowsybl `grid.xiidm`.

The environment used throughout is `case_14_matpower/`, whose `grid.m` is the IEEE 14
bus system exactly as `pandapower.networks.case14()` has it -- which is what makes the
comparison in `TestAgainstPandapower` an actual oracle rather than a tautology: the
same grid, read by two independent converters, solved by two independent solvers.
"""

import os
import tempfile
import unittest
import warnings

import numpy as np

import grid2op
from grid2op.Runner import Runner

from lightsim2grid import LightSimBackend
from lightsim2grid.network import init_from_matpower
from lightsim2grid.network.from_matpower._parse_matpower_source import load_matpower_data

try:
    from grid2op._create_test_suite import create_test_suite  # noqa: F401
    CAN_DO_TEST_SUITE = True
except ImportError:
    CAN_DO_TEST_SUITE = False

try:
    import pandapower as pp
    import pandapower.networks as pn
    CAN_DO_PANDAPOWER = True
except ImportError:
    CAN_DO_PANDAPOWER = False

try:
    import scipy.io  # noqa: F401
    CAN_DO_MAT = True
except ImportError:
    CAN_DO_MAT = False

try:
    import matpowercaseframes  # noqa: F401
    CAN_READ_DOT_M = True
except ImportError:
    # reading the ".m" a matpower case is distributed as is what this optional package
    # is for. Everything below reads `case_14_matpower/grid.m`, so without it there is
    # nothing to test rather than something failing
    CAN_READ_DOT_M = False


PATH_CASE_14_MATPOWER = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                     "case_14_matpower")
GRID_M = os.path.join(PATH_CASE_14_MATPOWER, "grid.m")

# how `case_14_matpower/grid.m` is laid out, read off the file itself (1-based
# matpower bus numbers, so one less than the grid2op substation ids below)
NB_SUB = 14
NB_LINE_ONLY = 17      # matpower "branch" rows with ratio == 0
NB_TRAFO = 3           # matpower "branch" rows with a ratio
NB_GEN = 5
NB_LOAD = 11
NB_SHUNT = 1


def _aux_prep_backend(backend, env_name):
    """What `grid2op.make` does to a backend, minus the environment around it."""
    type(backend)._clear_grid_dependant_class_attributes()
    backend.set_env_name(env_name)
    backend.load_grid(PATH_CASE_14_MATPOWER, "grid.m")
    backend.load_storage_data(PATH_CASE_14_MATPOWER)
    backend.load_redispacthing_data(PATH_CASE_14_MATPOWER)
    backend.assert_grid_correct()
    return backend


@unittest.skipIf(not CAN_READ_DOT_M, "matpowercaseframes is not installed")
class TestLoadGridMatpower(unittest.TestCase):
    """Reading `grid.m` gives a backend grid2op accepts, with the right shape."""

    def setUp(self) -> None:
        self.backend = _aux_prep_backend(LightSimBackend(loader_method="matpower"),
                                         type(self).__name__)

    def test_supported_grid_format(self):
        # this is what makes `grid2op.make(some_env)` look for a "grid.m" (then a
        # "grid.mat") instead of a "grid.json"
        assert self.backend.supported_grid_format == ("m", "mat")

    def test_shape(self):
        cls = type(self.backend)
        assert cls.n_sub == NB_SUB, f"wrong number of substations: {cls.n_sub} vs {NB_SUB}"
        assert cls.n_line == NB_LINE_ONLY + NB_TRAFO
        assert cls.n_gen == NB_GEN
        assert cls.n_load == NB_LOAD
        assert cls.n_shunt == NB_SHUNT
        assert cls.n_storage == 0, "matpower has no storage unit table"
        assert cls.shunts_data_available

    def test_one_substation_per_matpower_bus(self):
        # matpower has no notion of a busbar section within a bus, so a substation is
        # a bus, and the extra busbar sections grid2op asks for are added on top
        assert self.backend._LightSimBackend__nb_bus_before == NB_SUB
        assert self.backend.nb_bus_total == NB_SUB * type(self.backend).n_busbar_per_sub
        assert self.backend._LightSimBackend__nb_powerline == NB_LINE_ONLY

    def test_everything_starts_connected_on_busbar_1(self):
        assert (self.backend._LightSimBackend__init_topo_vect == 1).all()

    def test_names_are_unique(self):
        cls = type(self.backend)
        for attr_nm in ["name_sub", "name_load", "name_gen", "name_line", "name_shunt"]:
            names = getattr(cls, attr_nm)
            assert len(set(names)) == len(names), f"duplicated name in {attr_nm}: {names}"

    def test_names_are_grid2op_s_own_default_ones(self):
        """matpower names nothing, so grid2op's `Backend._fill_names_obj` makes the
        names -- the backend must not roll its own, or they would drift from what every
        other nameless grid gets."""
        cls = type(self.backend)
        reference = LightSimBackend(loader_method="matpower")
        reference.load_to_subid = cls.load_to_subid
        reference.gen_to_subid = cls.gen_to_subid
        reference.line_or_to_subid = cls.line_or_to_subid
        reference.line_ex_to_subid = cls.line_ex_to_subid
        reference.storage_to_subid = cls.storage_to_subid
        reference.shunt_to_subid = cls.shunt_to_subid
        reference.n_sub = cls.n_sub
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            reference._fill_names_obj()
        for attr_nm in ["name_sub", "name_load", "name_gen", "name_line", "name_shunt"]:
            np.testing.assert_array_equal(getattr(cls, attr_nm), getattr(reference, attr_nm),
                                          err_msg=f"{attr_nm} is not grid2op's default")
        # and the grid knows them too, so an LSGrid error names what grid2op names
        np.testing.assert_array_equal([el.name for el in self.backend._grid.get_loads()],
                                      cls.name_load)
        np.testing.assert_array_equal([el.name for el in self.backend._grid.get_generators()],
                                      cls.name_gen)

    def test_nominal_voltages(self):
        # the "baseKV" column of grid.m: 132 kV for buses 1-5, 11 kV for bus 8, 33 kV
        # for the rest -- and a substation's voltage is what its elements report
        cls = type(self.backend)
        vn_kv = np.array(self.backend._grid.get_bus_vn_kv())[:NB_SUB]
        expected = np.array([132., 132., 132., 132., 132., 33., 33., 11.,
                             33., 33., 33., 33., 33., 33.])
        np.testing.assert_allclose(vn_kv, expected)
        np.testing.assert_allclose(self.backend.prod_pu_to_kv, expected[cls.gen_to_subid])
        np.testing.assert_allclose(self.backend.load_pu_to_kv, expected[cls.load_to_subid])
        np.testing.assert_allclose(self.backend._sh_vnkv, expected[cls.shunt_to_subid])


@unittest.skipIf(not CAN_READ_DOT_M, "matpowercaseframes is not installed")
class TestSubidComeFromTheLoader(unittest.TestCase):
    """
    The substation each element sits in is a property of the file, so it is
    `init_from_matpower` (through the shared `init_from_powermodels` engine) that
    works it out and sets it on the `LSGrid`; the backend reads it back rather than
    deriving it a second time. These tests pin both halves of that.
    """

    def setUp(self) -> None:
        self.bus, self.gen, self.branch, _, _ = load_matpower_data(GRID_M)
        # matpower bus number -> lightsim2grid substation id
        self.mp_to_sub = {int(b): i for i, b in enumerate(self.bus[:, 0])}
        self.model = init_from_matpower(GRID_M)

    def _sub_of(self, matpower_bus_numbers):
        return np.array([self.mp_to_sub[int(b)] for b in matpower_bus_numbers])

    def test_loader_sets_them_on_the_lsgrid(self):
        model = self.model
        # generators: the "GEN_BUS" column of mpc.gen
        np.testing.assert_array_equal([el.sub_id for el in model.get_generators()],
                                      self._sub_of(self.gen[:, 0]))
        # loads: matpower has no load table, a load is a bus row with a non-zero PD/QD
        load_bus = self.bus[(self.bus[:, 2] != 0.) | (self.bus[:, 3] != 0.), 0]
        np.testing.assert_array_equal([el.sub_id for el in model.get_loads()],
                                      self._sub_of(load_bus))
        # shunts: a bus row with a non-zero GS/BS
        shunt_bus = self.bus[(self.bus[:, 4] != 0.) | (self.bus[:, 5] != 0.), 0]
        np.testing.assert_array_equal([el.sub_id for el in model.get_shunts()],
                                      self._sub_of(shunt_bus))
        # branches: ratio == 0 is a plain powerline, anything else a transformer, and
        # lightsim2grid keeps all the lines before all the transformers
        is_trafo = self.branch[:, 8] != 0.
        np.testing.assert_array_equal([el.sub1_id for el in model.get_lines()],
                                      self._sub_of(self.branch[~is_trafo, 0]))
        np.testing.assert_array_equal([el.sub2_id for el in model.get_lines()],
                                      self._sub_of(self.branch[~is_trafo, 1]))
        np.testing.assert_array_equal([el.sub1_id for el in model.get_trafos()],
                                      self._sub_of(self.branch[is_trafo, 0]))
        np.testing.assert_array_equal([el.sub2_id for el in model.get_trafos()],
                                      self._sub_of(self.branch[is_trafo, 1]))

    def test_backend_reads_them_back_unchanged(self):
        backend = _aux_prep_backend(LightSimBackend(loader_method="matpower"),
                                    type(self).__name__)
        cls = type(backend)
        model = self.model
        np.testing.assert_array_equal(cls.gen_to_subid, [el.sub_id for el in model.get_generators()])
        np.testing.assert_array_equal(cls.load_to_subid, [el.sub_id for el in model.get_loads()])
        np.testing.assert_array_equal(cls.shunt_to_subid, [el.sub_id for el in model.get_shunts()])
        np.testing.assert_array_equal(
            cls.line_or_to_subid,
            [el.sub1_id for el in model.get_lines()] + [el.sub1_id for el in model.get_trafos()])
        np.testing.assert_array_equal(
            cls.line_ex_to_subid,
            [el.sub2_id for el in model.get_lines()] + [el.sub2_id for el in model.get_trafos()])

    def test_an_out_of_service_element_still_knows_its_substation(self):
        """
        The reason the loader returns the bus it *built* each element on instead of
        reading `el.bus_id` back: a disconnected element reports `bus_id == -1`, while
        the substation it belongs to is a property of the grid, not of its status.
        """
        bus, gen, branch, _, baseMVA = load_matpower_data(GRID_M)
        gen = gen.copy()
        gen[2, 7] = 0  # GEN_STATUS of the bus-6 machine
        branch = branch.copy()
        branch[0, 10] = 0  # BR_STATUS of the branch 1 -- 2
        model = init_from_matpower({"bus": bus, "gen": gen, "branch": branch,
                                    "baseMVA": baseMVA})
        off_gen = [el for el in model.get_generators()][2]
        assert not off_gen.connected
        assert off_gen.bus_id == -1
        assert off_gen.sub_id == self.mp_to_sub[int(gen[2, 0])]
        off_line = [el for el in model.get_lines()][0]
        assert not off_line.connected1 and not off_line.connected2
        assert off_line.sub1_id == self.mp_to_sub[int(branch[0, 0])]
        assert off_line.sub2_id == self.mp_to_sub[int(branch[0, 1])]

    def test_a_topology_action_uses_them(self):
        """
        What the substation ids are actually *for*: `update_topo` turns "busbar 2 of
        my substation" into a global bus id with them, so a topology action landing on
        the right bus is the end to end proof they are right.
        """
        backend = _aux_prep_backend(LightSimBackend(loader_method="matpower"),
                                    type(self).__name__ + "_topo")
        cls = type(backend)
        action = type(backend)._complete_action_class()
        action.update({"set_bus": {"loads_id": [(3, 2)]}})
        bk_act = backend.my_bk_act_class()
        bk_act += action
        backend.apply_action(bk_act)
        load = [el for el in backend._grid.get_loads()][3]
        # busbar section 2 of substation `sub` is the global bus `sub + n_sub`
        assert load.bus_id == cls.load_to_subid[3] + NB_SUB, (
            f"load 3 landed on bus {load.bus_id} instead of "
            f"{cls.load_to_subid[3] + NB_SUB}")


@unittest.skipIf(not CAN_READ_DOT_M, "matpowercaseframes is not installed")
class TestRunpfMatpower(unittest.TestCase):
    def setUp(self) -> None:
        self.backend = _aux_prep_backend(LightSimBackend(loader_method="matpower"),
                                         type(self).__name__)

    def test_runpf_ac(self):
        conv, exc_ = self.backend.runpf(is_dc=False)
        assert conv, f"AC powerflow diverged: {exc_}"
        assert self.backend.can_output_theta

    def test_runpf_dc(self):
        conv, exc_ = self.backend.runpf(is_dc=True)
        assert conv, f"DC powerflow diverged: {exc_}"

    def test_same_as_the_bare_lsgrid(self):
        """The backend must not perturb the grid the loader built: the voltages it
        solves for are the ones a plain `init_from_matpower` + `ac_pf` gives."""
        conv, exc_ = self.backend.runpf(is_dc=False)
        assert conv, f"{exc_}"
        model = init_from_matpower(GRID_M)
        Vfinal = model.ac_pf(np.ones(NB_SUB, dtype=complex), 20, 1e-10)
        assert Vfinal.shape[0] == NB_SUB
        np.testing.assert_allclose(np.abs(self.backend.V[:NB_SUB]), np.abs(Vfinal), rtol=1e-8, atol=1e-8)


@unittest.skipIf(not CAN_DO_PANDAPOWER, "pandapower is not installed")
@unittest.skipIf(not CAN_READ_DOT_M, "matpowercaseframes is not installed")
class TestAgainstPandapower(unittest.TestCase):
    """
    `case_14_matpower/grid.m` holds the same IEEE 14 bus system as
    `pandapower.networks.case14()`. Read it through the matpower loader, solve it with
    lightsim2grid, and compare against pandapower solving its own copy: an independent
    converter and an independent solver, on the same grid.
    """

    def setUp(self) -> None:
        self.backend = _aux_prep_backend(LightSimBackend(loader_method="matpower"),
                                         type(self).__name__)
        conv, exc_ = self.backend.runpf(is_dc=False)
        assert conv, f"AC powerflow diverged: {exc_}"
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.net = pn.case14()
            pp.runpp(self.net, numba=False)

    def test_bus_voltages(self):
        # matpower bus `i` (1-based) is pandapower bus `i - 1` and lightsim2grid bus `i - 1`
        vm_ls = np.abs(self.backend.V[:NB_SUB])
        np.testing.assert_allclose(vm_ls, self.net.res_bus["vm_pu"].values, rtol=1e-6, atol=1e-6)

    def test_bus_angles(self):
        theta_ls = np.rad2deg(np.angle(self.backend.V[:NB_SUB]))
        np.testing.assert_allclose(theta_ls, self.net.res_bus["va_degree"].values, rtol=1e-5, atol=1e-5)

    def test_loads(self):
        load_p, load_q, _ = self.backend.loads_info()
        # both are ordered by (increasing) bus, so they line up element by element
        np.testing.assert_allclose(load_p, self.net.load["p_mw"].values, rtol=1e-6, atol=1e-6)
        np.testing.assert_allclose(load_q, self.net.load["q_mvar"].values, rtol=1e-6, atol=1e-6)

    def test_generation(self):
        prod_p, prod_q, _ = self.backend.generators_info()
        # pandapower splits the machines into 4 gens and one ext_grid (the slack); the
        # matpower case has the same 5 machines in the same order, the slack last
        pp_p = np.concatenate((self.net.res_gen["p_mw"].values, self.net.res_ext_grid["p_mw"].values))
        pp_q = np.concatenate((self.net.res_gen["q_mvar"].values, self.net.res_ext_grid["q_mvar"].values))
        np.testing.assert_allclose(prod_p, pp_p, rtol=1e-5, atol=1e-5)
        np.testing.assert_allclose(prod_q, pp_q, rtol=1e-5, atol=1e-5)

    def test_powerline_flows(self):
        p_or, q_or, _, _ = self.backend.lines_or_info()
        # the 15 matpower branches with r != 0 are pandapower's 15 lines, in order
        np.testing.assert_allclose(p_or[:15], self.net.res_line["p_from_mw"].values, rtol=1e-5, atol=1e-5)
        np.testing.assert_allclose(q_or[:15], self.net.res_line["q_from_mvar"].values, rtol=1e-5, atol=1e-5)

    def test_transformer_flows(self):
        p_or, q_or, _, _ = self.backend.lines_or_info()
        # matpower's 3 tap-changing transformers are pandapower's trafos 0, 1 and 2;
        # its 2 unity-ratio ones (branches "7 -- 8" and "7 -- 9", written with
        # matpower's "ratio == 0" plain-line sentinel) are pandapower's trafos 3 and 4
        # and are read as powerlines here, right after the 15 above
        np.testing.assert_allclose(p_or[NB_LINE_ONLY:], self.net.res_trafo["p_hv_mw"].values[:3],
                                   rtol=1e-5, atol=1e-5)
        np.testing.assert_allclose(p_or[15:NB_LINE_ONLY], self.net.res_trafo["p_hv_mw"].values[3:],
                                   rtol=1e-5, atol=1e-5)

    def test_shunt(self):
        sh_p, sh_q, _, _ = self.backend.shunt_info()
        np.testing.assert_allclose(sh_p, self.net.res_shunt["p_mw"].values, rtol=1e-5, atol=1e-5)
        np.testing.assert_allclose(sh_q, self.net.res_shunt["q_mvar"].values, rtol=1e-5, atol=1e-5)


@unittest.skipIf(not CAN_READ_DOT_M, "matpowercaseframes is not installed")
class TestLoaderKwargs(unittest.TestCase):
    def test_unknown_kwarg_raises(self):
        backend = LightSimBackend(loader_method="matpower",
                                  loader_kwargs={"i_do_not_exist": True})
        with self.assertRaises(RuntimeError):
            backend.load_grid(PATH_CASE_14_MATPOWER, "grid.m")

    def test_pypowsybl_only_kwarg_raises(self):
        # the matpower loader has its own (much smaller) set of accepted keys
        backend = LightSimBackend(loader_method="matpower",
                                  loader_kwargs={"use_buses_for_sub": True})
        with self.assertRaises(RuntimeError):
            backend.load_grid(PATH_CASE_14_MATPOWER, "grid.m")

    def test_grid_kwarg_skips_the_file(self):
        """`loader_kwargs={"grid": <an already parsed case>}` uses that case and
        ignores the path, exactly like the pypowsybl loader's own "grid" kwarg."""
        bus, gen, branch, _, baseMVA = load_matpower_data(GRID_M)
        mpc = {"bus": bus, "gen": gen, "branch": branch, "baseMVA": baseMVA}
        backend = LightSimBackend(loader_method="matpower", loader_kwargs={"grid": mpc})
        type(backend)._clear_grid_dependant_class_attributes()
        backend.set_env_name(type(self).__name__ + "_grid_kwarg")
        # a path that does not exist: it must not even be looked at
        backend.load_grid(PATH_CASE_14_MATPOWER, "i_am_not_a_file.m")
        backend.load_storage_data(PATH_CASE_14_MATPOWER)
        backend.load_redispacthing_data(PATH_CASE_14_MATPOWER)
        backend.assert_grid_correct()
        assert type(backend).n_sub == NB_SUB
        conv, exc_ = backend.runpf()
        assert conv, f"{exc_}"

    def test_n_busbar_per_sub(self):
        backend = LightSimBackend(loader_method="matpower",
                                  loader_kwargs={"n_busbar_per_sub": 3})
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            _aux_prep_backend(backend, type(self).__name__ + "_nbb3")
        assert type(backend).n_busbar_per_sub == 3
        assert backend.nb_bus_total == 3 * NB_SUB

    def test_double_bus_per_sub_is_not_a_matpower_kwarg(self):
        # it predates grid2op supporting any number of busbars per substation; only
        # the pypowsybl loader keeps it, for backward compatibility
        backend = LightSimBackend(loader_method="matpower",
                                  loader_kwargs={"double_bus_per_sub": True})
        with self.assertRaises(RuntimeError):
            backend.load_grid(PATH_CASE_14_MATPOWER, "grid.m")


@unittest.skipIf(not CAN_DO_MAT or not CAN_READ_DOT_M,
                 "scipy or matpowercaseframes is not installed")
class TestMatFile(unittest.TestCase):
    """An environment can just as well ship the ".mat" matpower saves."""

    def test_load_a_mat_grid(self):
        bus, gen, branch, _, baseMVA = load_matpower_data(GRID_M)
        with tempfile.TemporaryDirectory() as tmp_dir:
            scipy.io.savemat(os.path.join(tmp_dir, "grid.mat"),
                             {"mpc": {"bus": bus, "gen": gen, "branch": branch,
                                      "baseMVA": baseMVA, "version": "2"}})
            backend = LightSimBackend(loader_method="matpower")
            type(backend)._clear_grid_dependant_class_attributes()
            backend.set_env_name(type(self).__name__)
            backend.load_grid(tmp_dir, "grid.mat")
            backend.load_storage_data(tmp_dir)
            backend.load_redispacthing_data(tmp_dir)
            backend.assert_grid_correct()
            assert type(backend).n_sub == NB_SUB
            assert type(backend).n_line == NB_LINE_ONLY + NB_TRAFO
            conv, exc_ = backend.runpf()
            assert conv, f"{exc_}"


@unittest.skipIf(not CAN_READ_DOT_M, "matpowercaseframes is not installed")
class TestEnvMatpower(unittest.TestCase):
    """The point of all this: a real grid2op environment on a `grid.m` file."""

    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make(PATH_CASE_14_MATPOWER,
                                    backend=LightSimBackend(loader_method="matpower"),
                                    test=True,
                                    _add_to_name=type(self).__name__)
        super().setUp()

    def tearDown(self) -> None:
        self.env.close()
        return super().tearDown()

    def test_can_make_and_reset(self):
        obs = self.env.reset()
        assert obs.n_sub == NB_SUB
        assert obs.n_line == NB_LINE_ONLY + NB_TRAFO
        # the thermal limits of the environment's config.py are the base case flows
        # with a 50% margin: not "infinite", and not violated
        assert 0.5 < obs.rho.max() < 1.0, f"unexpected max rho: {obs.rho.max()}"

    def test_can_step(self):
        self.env.reset()
        obs, reward, done, info = self.env.step(self.env.action_space())
        assert not done, f"do-nothing ended the episode: {info['exception']}"
        assert not info["exception"]

    def test_topology_action(self):
        """
        Move every element of substation 0 (the origin of powerlines 0 and 1, and the
        reference-bus machine, generator 4) onto busbar 2. That is a pure relabelling
        of the substation, so the powerflow must give exactly the same answer -- which
        makes it a check of the substation ids the loader wrote, not just of the plumbing.
        """
        obs_before = self.env.reset()
        act = self.env.action_space({"set_bus": {"lines_or_id": [(0, 2), (1, 2)],
                                                 "generators_id": [(4, 2)]}})
        obs, reward, done, info = self.env.step(act)
        assert not done, f"the topology action ended the episode: {info['exception']}"
        assert not info["exception"]
        cls = type(self.env.backend)
        assert obs.topo_vect[cls.line_or_pos_topo_vect[0]] == 2
        assert obs.topo_vect[cls.line_or_pos_topo_vect[1]] == 2
        assert obs.topo_vect[cls.gen_pos_topo_vect[4]] == 2
        np.testing.assert_allclose(obs.p_or, obs_before.p_or, rtol=1e-6, atol=1e-6)
        np.testing.assert_allclose(obs.gen_p, obs_before.gen_p, rtol=1e-6, atol=1e-6)

    def test_disconnect_a_line(self):
        self.env.reset()
        obs, reward, done, info = self.env.step(
            self.env.action_space({"set_line_status": [(3, -1)]}))
        assert not done, f"the disconnection ended the episode: {info['exception']}"
        assert not obs.line_status[3]

    def test_copy(self):
        obs = self.env.reset()
        env_cpy = self.env.copy()
        try:
            obs_cpy = env_cpy.reset()
            assert self.env.backend.supported_grid_format == ("m", "mat")
            assert env_cpy.backend.supported_grid_format == ("m", "mat")
            assert env_cpy.backend._loader_method == "matpower"
            np.testing.assert_allclose(obs.rho, obs_cpy.rho)
        finally:
            env_cpy.close()

    def test_runner(self):
        self.env.reset()
        env_cpy = self.env.copy()
        try:
            runner = Runner(**self.env.get_params_for_runner())
            runner_cpy = Runner(**env_cpy.get_params_for_runner())
            res = runner.run(nb_episode=1, max_iter=5)
            res_cpy = runner_cpy.run(nb_episode=1, max_iter=5)
            assert len(res) == 1
            # `ChangeNothing` stops one step short of `max_iter`; what matters here is
            # that nothing ended the episode early (a game over would show up as a
            # much smaller step count)
            assert res[0][3] >= 4, f"the episode stopped early: {res[0]}"
            for el, el_cpy in zip(res[0], res_cpy[0]):
                assert el == el_cpy, f"{el} vs {el_cpy}"
        finally:
            env_cpy.close()


if CAN_DO_TEST_SUITE and CAN_READ_DOT_M:
    from grid2op.tests.aaa_test_backend_interface import AAATestBackendAPI

    class TestBackendAPI_MatpowerBk(AAATestBackendAPI, unittest.TestCase):
        """grid2op's own "does this backend respect the Backend API" test suite."""

        def get_path(self):
            return PATH_CASE_14_MATPOWER

        def get_casefile(self):
            return "grid.m"

        def make_backend(self, detailed_infos_for_cascading_failures=False):
            return LightSimBackend(
                loader_method="matpower",
                detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)

        def test_01load_grid(self):
            """
            Skipped on purpose. This one test of the suite asserts the *shape* of
            grid2op's own `educ_case14_storage` ("This test will NOT pass if the grid
            is not the educ_case14_storage file", says its docstring): 6 generators and
            2 storage units. `case_14_matpower/grid.m` is the IEEE 14 bus case, which
            has 5 machines and -- matpower having no storage table at all -- no storage
            unit. `TestLoadGridMatpower` above checks the same thing against the shape
            this grid actually has. Every other test of the suite runs.
            """
            self.skipTest("the grid is not grid2op's educ_case14_storage, which this "
                          "particular test of the suite hard-codes the shape of")


if __name__ == "__main__":
    unittest.main()
