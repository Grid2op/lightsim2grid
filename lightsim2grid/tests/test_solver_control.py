# TODO: test that "if I do something, then with or without the speed optim, I have the same thing"

# use backend._grid.tell_solver_need_reset() to force the reset of the solver
# and by "something" do :
# - disconnect line X
# - disconnect trafo X
# - reconnect line X
# - reconnect trafo X
# - change load X
# - change gen X
# - change shunt X
# - change storage X
# - change slack X
# - change slack weight X
# - turnoff_gen_pv X
# - when it diverges (and then I can make it un converge normally) X
# - test change bus to -1 and deactivate element has the same impact (irrelevant !)
# - test when set_bus to 2

# TODO and do that for all solver type: NR X, NRSingleSlack, GaussSeidel, GaussSeidelSynch, FDPF_XB and FDPFBX
# and for all solver check everything that can be check: final tolerance, number of iteration, etc. etc.

import unittest
import warnings
import numpy as np
import grid2op
from grid2op.Action import CompleteAction
from grid2op.Parameters import Parameters

from lightsim2grid import LightSimBackend
from lightsim2grid.solver import SolverType


class TestSolverControl(unittest.TestCase):
    def _aux_setup_grid(self):
        self.need_dc = True  # is it worth it to run DC powerflow ?
        self.can_dist_slack = True
        self.gridmodel.change_solver(SolverType.SparseLU)
        self.gridmodel.change_solver(SolverType.DC)
        
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("educ_case14_storage",
                                    test=True,
                                    action_class=CompleteAction,
                                    backend=LightSimBackend())
        self.gridmodel = self.env.backend._grid
        self.iter = 10
        self._aux_setup_grid()
        self.v_init = 0.0 * self.env.backend.V + 1.04  # just to have a vector with the right dimension
        self.tol_solver = 1e-8  # solver
        self.tol_equal = 1e-10  #  for comparing with and without the "smarter solver" things, and make sure everything is really equal!
    
    def test_update_topo_ac(self, runpf_fun="_run_ac_pf"):
        """test when I disconnect a line alone at one end: it changes the size of the ybus / sbus vector AC"""
        LINE_ID = 2
        dim_topo = type(self.env).dim_topo
        mask_changed = np.zeros(dim_topo, dtype=bool)
        mask_val = np.zeros(dim_topo, dtype=np.int32)
        mask_changed[type(self.env).line_ex_pos_topo_vect[LINE_ID]] = True
        mask_val[type(self.env).line_ex_pos_topo_vect[LINE_ID]] = 2
        self.gridmodel.update_topo(mask_changed, mask_val)
        V = getattr(self, runpf_fun)(gridmodel=self.gridmodel)
        assert len(V), "it should not have diverged here"
        self.gridmodel.unset_changes()
        
        mask_changed = np.zeros(dim_topo, dtype=bool)
        mask_val = np.zeros(dim_topo, dtype=np.int32)
        mask_changed[type(self.env).line_or_pos_topo_vect[LINE_ID]] = True
        mask_val[type(self.env).line_or_pos_topo_vect[LINE_ID]] = -1
        # mask_changed[type(self.env).line_ex_pos_topo_vect[LINE_ID]] = True
        # mask_val[type(self.env).line_ex_pos_topo_vect[LINE_ID]] = -1
        self.gridmodel.update_topo(mask_changed, mask_val)
        solver = self.gridmodel.get_dc_solver() if runpf_fun == "_run_dc_pf" else self.gridmodel.get_solver() 
        V1 = getattr(self, runpf_fun)(gridmodel=self.gridmodel)
        assert len(V1), f"it should not have diverged here. Error : {solver.get_error()}"
        
        self.gridmodel.tell_solver_need_reset()
        V2 = getattr(self, runpf_fun)(gridmodel=self.gridmodel)
        assert len(V2), f"it should not have diverged here. Error : {solver.get_error()}"
        assert np.allclose(V1, V2, rtol=self.tol_equal, atol=self.tol_equal)
    
    def test_update_topo_dc(self):
        """test when I disconnect a line alone at one end: it changes the size of the ybus / sbus vector AC"""
        if not self.need_dc:
            self.skipTest("Useless to run DC")
        self.test_update_topo_ac("_run_dc_pf")   
        
    def test_pf_run_dc(self):
        """test I have the same results if nothing is done with and without restarting from scratch when running dc powerflow"""
        if not self.need_dc:
            self.skipTest("Useless to run DC")
        Vdc_init = self.gridmodel.dc_pf(self.v_init, self.iter, self.tol_solver)
        assert len(Vdc_init), f"error: gridmodel should converge in DC"
        self.gridmodel.unset_changes()
        Vdc_init2 = self.gridmodel.dc_pf(self.v_init, self.iter, self.tol_solver)
        assert len(Vdc_init2), f"error: gridmodel should converge in DC"
        self.gridmodel.tell_solver_need_reset()
        Vdc_init3 = self.gridmodel.dc_pf(self.v_init, self.iter, self.tol_solver)
        assert len(Vdc_init3), f"error: gridmodel should converge in DC"
        assert np.allclose(Vdc_init, Vdc_init2, rtol=self.tol_equal, atol=self.tol_equal)
        assert np.allclose(Vdc_init2, Vdc_init3, rtol=self.tol_equal, atol=self.tol_equal)
        
    def test_pf_run_ac(self):   
        """test I have the same results if nothing is done with and without restarting from scratch when running ac powerflow"""
        Vac_init = self.gridmodel.ac_pf(self.v_init, self.iter, self.tol_solver)
        assert len(Vac_init), f"error: gridmodel should converge in AC"
        self.gridmodel.unset_changes()
        Vac_init2 = self.gridmodel.ac_pf(self.v_init, self.iter, self.tol_solver)
        assert len(Vac_init2), f"error: gridmodel should converge in AC"
        self.gridmodel.tell_solver_need_reset()
        Vac_init3 = self.gridmodel.ac_pf(self.v_init, self.iter, self.tol_solver)
        assert len(Vac_init3), f"error: gridmodel should converge in AC"
        assert np.allclose(Vac_init, Vac_init2, rtol=self.tol_equal, atol=self.tol_equal)
        assert np.allclose(Vac_init2, Vac_init3, rtol=self.tol_equal, atol=self.tol_equal)
    
    def _disco_line_action(self, gridmodel, el_id=0, el_val=0.):
        gridmodel.deactivate_powerline(el_id)
        
    def _reco_line_action(self, gridmodel, el_id=0, el_val=0.):
        gridmodel.reactivate_powerline(el_id)
    
    def _disco_trafo_action(self, gridmodel, el_id=0, el_val=0.):
        gridmodel.deactivate_trafo(el_id)
        
    def _reco_trafo_action(self, gridmodel, el_id=0, el_val=0.):
        gridmodel.reactivate_trafo(el_id)
    
    def _run_ac_pf(self, gridmodel):
        return gridmodel.ac_pf(self.v_init, self.iter, self.tol_solver)
    
    def _run_dc_pf(self, gridmodel):
        return gridmodel.dc_pf(self.v_init, self.iter, self.tol_solver)
    
    def aux_do_undo_ac(self,
                       runpf_fun="_run_ac_pf",
                       funname_do="_disco_line_action",
                       funname_undo="_reco_line_action",
                       el_id=0,
                       el_val=0.,
                       to_add_remove=0.,
                       expected_diff=0.1
                       ):
        pf_mode = "AC" if runpf_fun=="_run_ac_pf" else "DC"
                
        # test "do the action"
        V_init = getattr(self, runpf_fun)(gridmodel=self.gridmodel)
        assert len(V_init), f"error for el_id={el_id}: gridmodel should converge in {pf_mode}"
        getattr(self, funname_do)(gridmodel=self.gridmodel, el_id=el_id, el_val=el_val + to_add_remove)
        V_disc = getattr(self, runpf_fun)(gridmodel=self.gridmodel)
        if len(V_disc) > 0:
            # powerflow converges, all should converge
            assert (np.abs(V_init - V_disc) >= expected_diff).any(), f"error for el_id={el_id}: at least one bus should have changed its result voltage in {pf_mode}: max {np.abs(V_init - V_disc).max():.2e}"
            self.gridmodel.unset_changes()
            V_disc1 = getattr(self, runpf_fun)(gridmodel=self.gridmodel)
            assert len(V_disc1), f"error for el_id={el_id}: gridmodel should converge in {pf_mode}"
            assert np.allclose(V_disc, V_disc1, rtol=self.tol_equal, atol=self.tol_equal)
            self.gridmodel.tell_solver_need_reset()
            V_disc2 = getattr(self, runpf_fun)(gridmodel=self.gridmodel)
            assert len(V_disc2), f"error for el_id={el_id}: gridmodel should converge in {pf_mode}"
            assert np.allclose(V_disc2, V_disc1, rtol=self.tol_equal, atol=self.tol_equal)
            assert np.allclose(V_disc2, V_disc, rtol=self.tol_equal, atol=self.tol_equal)
        else:
            #powerflow diverges
            self.gridmodel.unset_changes()
            V_disc1 = getattr(self, runpf_fun)(gridmodel=self.gridmodel)
            assert len(V_disc1) == 0, f"error for el_id={el_id}: powerflow should diverge as it did initially in {pf_mode}"
            self.gridmodel.tell_solver_need_reset()
            V_disc2 = getattr(self, runpf_fun)(gridmodel=self.gridmodel)
            assert len(V_disc1) == 0, f"error for el_id={el_id}: powerflow should diverge as it did initially in {pf_mode}"
            
        # test "undo the action"
        self.gridmodel.unset_changes()
        getattr(self, funname_undo)(gridmodel=self.gridmodel, el_id=el_id, el_val=el_val)
        V_reco = getattr(self, runpf_fun)(gridmodel=self.gridmodel)
        assert len(V_reco), f"error for el_id={el_id}: gridmodel should converge in {pf_mode}"
        assert np.allclose(V_reco, V_init, rtol=self.tol_equal, atol=self.tol_equal), f"error for el_id={el_id}: do an action and then undo it should not have any impact in {pf_mode}: max {np.abs(V_init - V_reco).max():.2e}"
        self.gridmodel.unset_changes()
        V_reco1 = getattr(self, runpf_fun)(gridmodel=self.gridmodel)
        assert len(V_reco1), f"error for el_id={el_id}: gridmodel should converge in {pf_mode}"
        assert np.allclose(V_reco1, V_reco, rtol=self.tol_equal, atol=self.tol_equal)
        self.gridmodel.tell_solver_need_reset()
        V_reco2 = getattr(self, runpf_fun)(gridmodel=self.gridmodel)
        assert len(V_reco2), f"error for el_id={el_id}: gridmodel should converge in {pf_mode}"
        assert np.allclose(V_reco2, V_reco1, rtol=self.tol_equal, atol=self.tol_equal)
        assert np.allclose(V_reco2, V_reco, rtol=self.tol_equal, atol=self.tol_equal)
            
    def test_disco_reco_line_ac(self, runpf_fun="_run_ac_pf"):   
        """test I have the same results if I disconnect a line with and 
        without restarting from scratch when running ac powerflow"""
        for el_id in range(len(self.gridmodel.get_lines())):
            self.gridmodel.tell_solver_need_reset()
            expected_diff = 3e-2
            if runpf_fun=="_run_ac_pf":
                if el_id == 4:
                    expected_diff = 1e-2
                elif el_id == 10:
                    expected_diff = 1e-3
                elif el_id == 12:
                    expected_diff = 1e-2
                elif el_id == 13:
                    expected_diff = 3e-3
            elif runpf_fun=="_run_dc_pf":
                if el_id == 4:
                    expected_diff = 1e-2
                elif el_id == 8:
                    expected_diff = 1e-2
                elif el_id == 10:
                    expected_diff = 1e-2
                elif el_id == 12:
                    expected_diff = 1e-2
                elif el_id == 13:
                    expected_diff = 3e-3
                elif el_id == 14:
                    expected_diff = 1e-2
            self.aux_do_undo_ac(funname_do="_disco_line_action",
                                funname_undo="_reco_line_action",
                                runpf_fun=runpf_fun,
                                el_id=el_id,
                                expected_diff=expected_diff
                                )
            
    def test_disco_reco_line_dc(self):   
        """test I have the same results if I disconnect a line with and 
        without restarting from scratch when running DC powerflow"""
        if not self.need_dc:
            self.skipTest("Useless to run DC")
        self.test_disco_reco_line_ac(runpf_fun="_run_dc_pf")

    def test_disco_reco_trafo_ac(self, runpf_fun="_run_ac_pf"):   
        """test I have the same results if I disconnect a trafo with and 
        without restarting from scratch when running ac powerflow"""
        for el_id in range(len(self.gridmodel.get_trafos())):
            self.gridmodel.tell_solver_need_reset()
            expected_diff = 3e-2
            if runpf_fun=="_run_ac_pf":
                if el_id == 1:
                    expected_diff = 1e-2
            self.aux_do_undo_ac(funname_do="_disco_trafo_action",
                                funname_undo="_reco_trafo_action",
                                runpf_fun=runpf_fun,
                                el_id=el_id,
                                expected_diff=expected_diff
                                )
            
    def test_disco_reco_trafo_dc(self):   
        """test I have the same results if I disconnect a trafo with and 
        without restarting from scratch when running DC powerflow"""
        if not self.need_dc:
            self.skipTest("Useless to run DC")
        self.test_disco_reco_trafo_ac(runpf_fun="_run_dc_pf")

    def _change_load_p_action(self, gridmodel, el_id, el_val):
        gridmodel.change_p_load(el_id, el_val)
        
    def test_change_load_p_ac(self, runpf_fun="_run_ac_pf"):   
        """test I have the same results if I change the load p with and 
        without restarting from scratch when running ac powerflow"""
        for load in self.gridmodel.get_loads():
            self.gridmodel.tell_solver_need_reset()
            expected_diff = 1e-3
            to_add_remove = 5.
            self.aux_do_undo_ac(funname_do="_change_load_p_action",
                                funname_undo="_change_load_p_action",
                                runpf_fun=runpf_fun,
                                el_id=load.id,
                                expected_diff=expected_diff,
                                el_val=load.target_p_mw,
                                to_add_remove=to_add_remove,
                                )
            
    def test_change_load_p_dc(self):   
        """test I have the same results if I change the load p with and 
        without restarting from scratch when running dc powerflow"""
        if not self.need_dc:
            self.skipTest("Useless to run DC")
        self.test_change_load_p_ac(runpf_fun="_run_dc_pf")

    def _change_load_q_action(self, gridmodel, el_id, el_val):
        gridmodel.change_q_load(el_id, el_val)
        
    def _change_gen_p_action(self, gridmodel, el_id, el_val):
        gridmodel.change_p_gen(el_id, el_val)
        
    def test_change_load_q_ac(self, runpf_fun="_run_ac_pf"):   
        """test I have the same results if I change the load q with and 
        without restarting from scratch when running ac powerflow
        
        NB: this test has no sense in DC
        """
        gen_bus = [el.bus_id for el in self.gridmodel.get_generators()]
        for load in self.gridmodel.get_loads():
            if load.bus_id in gen_bus:
                # nothing will change if there is a gen connected to the
                # same bus as the load
                continue
            self.gridmodel.tell_solver_need_reset()
            expected_diff = 1e-3
            to_add_remove = 5.
            self.aux_do_undo_ac(funname_do="_change_load_q_action",
                                funname_undo="_change_load_q_action",
                                runpf_fun=runpf_fun,
                                el_id=load.id,
                                expected_diff=expected_diff,
                                el_val=load.target_q_mvar,
                                to_add_remove=to_add_remove,
                                )
            
    def test_change_gen_p_ac(self, runpf_fun="_run_ac_pf"):   
        """test I have the same results if I change the gen p with and 
        without restarting from scratch when running ac powerflow"""
        for gen in self.gridmodel.get_generators():
            if gen.is_slack:
                # nothing to do for the slack bus...
                continue
            self.gridmodel.tell_solver_need_reset()
            expected_diff = 1e-3
            to_add_remove = 5.
            self.aux_do_undo_ac(funname_do="_change_gen_p_action",
                                funname_undo="_change_gen_p_action",
                                runpf_fun=runpf_fun,
                                el_id=gen.id,
                                expected_diff=expected_diff,
                                el_val=gen.target_p_mw,
                                to_add_remove=to_add_remove,
                                )
            
    def test_change_gen_p_dc(self):   
        """test I have the same results if I change the gen p with and 
        without restarting from scratch when running dc powerflow"""
        if not self.need_dc:
            self.skipTest("Useless to run DC")
        self.test_change_gen_p_ac(runpf_fun="_run_dc_pf")
        
    def _change_gen_v_action(self, gridmodel, el_id, el_val):
        gridmodel.change_v_gen(el_id, el_val)
        
    def test_change_gen_v_ac(self, runpf_fun="_run_ac_pf"):   
        """test I have the same results if I change the gen v with and 
        without restarting from scratch when running ac powerflow
        
        NB: this test has no sense in DC
        """
        gen_bus = [el.bus_id for el in self.gridmodel.get_generators()]
        vn_kv = self.gridmodel.get_bus_vn_kv()
        for gen in self.gridmodel.get_generators():
            if gen.bus_id in gen_bus:
                # nothing will change if there is another gen connected to the
                # same bus as the gen (error if everything is coded normally
                # which it might not)
                continue
            
            self.gridmodel.tell_solver_need_reset()
            expected_diff = 1e-3
            to_add_remove = 0.1 * gen.target_vm_pu * vn_kv[gen.bus_id]
            self.aux_do_undo_ac(funname_do="_change_gen_v_action",
                                funname_undo="_change_gen_v_action",
                                runpf_fun=runpf_fun,
                                el_id=gen.id,
                                expected_diff=expected_diff,
                                el_val=gen.target_vm_pu * vn_kv[gen.bus_id],
                                to_add_remove=to_add_remove,
                                )

    def _change_shunt_p_action(self, gridmodel, el_id, el_val):
        gridmodel.change_p_shunt(el_id, el_val)
        
    def test_change_shunt_p_ac(self, runpf_fun="_run_ac_pf"):   
        """test I have the same results if I change the shunt p with and 
        without restarting from scratch when running ac powerflow"""
        for shunt in self.gridmodel.get_shunts():
            self.gridmodel.tell_solver_need_reset()
            expected_diff = 1e-3
            to_add_remove = 5.
            self.aux_do_undo_ac(funname_do="_change_shunt_p_action",
                                funname_undo="_change_shunt_p_action",
                                runpf_fun=runpf_fun,
                                el_id=shunt.id,
                                expected_diff=expected_diff,
                                el_val=shunt.target_p_mw,
                                to_add_remove=to_add_remove,
                                )
            
    def test_change_shunt_p_dc(self):   
        """test I have the same results if I change the shunt p with and 
        without restarting from scratch when running dc powerflow"""
        if not self.need_dc:
            self.skipTest("Useless to run DC")
        self.test_change_shunt_p_ac(runpf_fun="_run_dc_pf")

    def _change_shunt_q_action(self, gridmodel, el_id, el_val):
        gridmodel.change_q_shunt(el_id, el_val)
        
    def test_change_shunt_q_ac(self, runpf_fun="_run_ac_pf"):   
        """test I have the same results if I change the shunt q with and 
        without restarting from scratch when running ac powerflow
        
        NB dc is not needed here (and does not make sense)"""
        for shunt in self.gridmodel.get_shunts():
            self.gridmodel.tell_solver_need_reset()
            expected_diff = 1e-3
            to_add_remove = 5.
            self.aux_do_undo_ac(funname_do="_change_shunt_q_action",
                                funname_undo="_change_shunt_q_action",
                                runpf_fun=runpf_fun,
                                el_id=shunt.id,
                                expected_diff=expected_diff,
                                el_val=shunt.target_q_mvar,
                                to_add_remove=to_add_remove,
                                )

    def _change_storage_p_action(self, gridmodel, el_id, el_val):
        gridmodel.change_p_storage(el_id, el_val)
        
    def test_change_storage_p_ac(self, runpf_fun="_run_ac_pf"):   
        """test I have the same results if I change the storage p with and 
        without restarting from scratch when running ac powerflow"""
        for storage in self.gridmodel.get_storages():
            self.gridmodel.tell_solver_need_reset()
            expected_diff = 1e-3
            to_add_remove = 5.
            self.aux_do_undo_ac(funname_do="_change_storage_p_action",
                                funname_undo="_change_storage_p_action",
                                runpf_fun=runpf_fun,
                                el_id=storage.id,
                                expected_diff=expected_diff,
                                el_val=storage.target_p_mw,
                                to_add_remove=to_add_remove,
                                )
            
    def test_change_storage_p_dc(self):   
        """test I have the same results if I change the storage p with and 
        without restarting from scratch when running dc powerflow"""
        if not self.need_dc:
            self.skipTest("Useless to run DC")
        self.test_change_storage_p_ac(runpf_fun="_run_dc_pf")

    def _change_storage_q_action(self, gridmodel, el_id, el_val):
        gridmodel.change_q_storage(el_id, el_val)
        
    def test_change_storage_q_ac(self, runpf_fun="_run_ac_pf"):   
        """test I have the same results if I change the storage q with and 
        without restarting from scratch when running ac powerflow"""
        gen_bus = [el.bus_id for el in self.gridmodel.get_generators()]
        for storage in self.gridmodel.get_storages():
            if storage.bus_id in gen_bus:
                # nothing will change if there is a gen connected to the
                # same bus as the load
                continue
            self.gridmodel.tell_solver_need_reset()
            expected_diff = 1e-3
            to_add_remove = 5.
            self.aux_do_undo_ac(funname_do="_change_storage_q_action",
                                funname_undo="_change_storage_q_action",
                                runpf_fun=runpf_fun,
                                el_id=storage.id,
                                expected_diff=expected_diff,
                                el_val=storage.target_q_mvar,
                                to_add_remove=to_add_remove,
                                )

    def _change_unique_slack_id(self, gridmodel, el_id, el_val):
        # el_id : new slack
        # el_val : old slack
        gridmodel.remove_gen_slackbus(el_val)  # remove old slack
        gridmodel.add_gen_slackbus(el_id, 1.0)  # add new slack

    def _unchange_unique_slack_id(self, gridmodel, el_id, el_val):
        # el_id : new slack
        # el_val : old slack
        gridmodel.remove_gen_slackbus(el_id)  # remove new slack
        gridmodel.add_gen_slackbus(el_val, 1.0)  # add back old slack

    def test_change_gen_slack_unique_ac(self, runpf_fun="_run_ac_pf"):   
        """test I have the same results if I change the (for this test) unique slack bus
        This is done in AC"""
        gen_is_slack = [el.is_slack for el in self.gridmodel.get_generators()]
        gen_id_slack_init = np.where(gen_is_slack)[0][0]
        
        for gen in self.gridmodel.get_generators():
            if gen_is_slack[gen.id]:
                # generator was already slack, don't expect any change
                continue
            self.gridmodel.tell_solver_need_reset()
            expected_diff = 1e-5  # change slack impact is low
            self.aux_do_undo_ac(funname_do="_change_unique_slack_id",
                                funname_undo="_unchange_unique_slack_id",
                                runpf_fun=runpf_fun,
                                el_id=gen.id,  # new slack
                                expected_diff=expected_diff,
                                el_val=gen_id_slack_init,  # old slack
                                to_add_remove=0,
                                )
            
    def test_change_gen_slack_unique_dc(self):   
        """test I have the same results if I change the (for this test) unique slack bus
        This is done in DC"""
        if not self.need_dc:
            self.skipTest("Useless to run DC")
        self.test_change_gen_slack_unique_ac(runpf_fun="_run_dc_pf")

    def _change_dist_slack_id(self, gridmodel, el_id, el_val):
        # el_id : new slack
        # el_val : old slack
        gridmodel.add_gen_slackbus(el_val, 0.5)  # set all slack a weight of 0.5
        gridmodel.add_gen_slackbus(el_id, 0.5)  # add new slack with a weight of 0.5

    def _unchange_dist_slack_id(self, gridmodel, el_id, el_val):
        # el_id : new slack
        # el_val : old slack
        gridmodel.remove_gen_slackbus(el_id)  # remove new slack
        gridmodel.add_gen_slackbus(el_val, 1.0)  # add back old slack with a weight of 1. (as originally)
    
    def test_distslack_weight_ac(self):
        """test I have the same results if the slack weights (dist mode) are changed
        
        NB: not done in DC because as of now DC does not support dist slack
        """
        if not self.can_dist_slack :
            self.skipTest("This solver does not support distributed slack")
        gen_is_slack = [el.is_slack for el in self.gridmodel.get_generators()]
        gen_id_slack_init = np.where(gen_is_slack)[0][0]
        runpf_fun = "_run_ac_pf"
        for gen in self.gridmodel.get_generators():
            if gen_is_slack[gen.id]:
                # generator was already slack, don't expect any change
                continue
            self.gridmodel.tell_solver_need_reset()
            expected_diff = 1e-5  # change slack impact is low
            self.aux_do_undo_ac(funname_do="_change_dist_slack_id",
                                funname_undo="_unchange_dist_slack_id",
                                runpf_fun=runpf_fun,
                                el_id=gen.id,  # new added slack
                                expected_diff=expected_diff,
                                el_val=gen_id_slack_init,  # old slack
                                to_add_remove=0,
                                )

    def _change_turnedoff_pv(self, gridmodel, el_id, el_val):
        gridmodel.turnedoff_no_pv()
        
    def _unchange_turnedoff_pv(self, gridmodel, el_id, el_val):
        gridmodel.turnedoff_pv()
    
    def test_turnedoff_pv_ac(self):
        """test the `turnedoff_pv` functionality
        """
        runpf_fun = "_run_ac_pf"
        self.gridmodel.tell_solver_need_reset()
        expected_diff = 1e-2  # change slack impact is low
        self.aux_do_undo_ac(funname_do="_change_turnedoff_pv",
                            funname_undo="_unchange_turnedoff_pv",
                            runpf_fun=runpf_fun,
                            el_id=0,
                            expected_diff=expected_diff,
                            el_val=0,
                            to_add_remove=0,
                            )

    def _change_for_divergence_sbus(self, gridmodel, el_id, el_val):
        # el_id=load_p_init,  # initial loads
        # el_val=load_p_div,  # final loads
        [gridmodel.change_p_load(l_id, val) for l_id, val in enumerate(el_val)]
        
    def _unchange_for_divergence_sbus(self, gridmodel, el_id, el_val):
        # el_id=load_p_init,  # initial loads
        # el_val=load_p_div,  # final loads
        [gridmodel.change_p_load(l_id, val) for l_id, val in enumerate(el_id)]
            
    def test_divergence_sbus_ac(self):
        """test I can make the grid diverge and converge again using sbus (ac mode only)
        
        It never diverge in DC because of Sbus"""
        runpf_fun = "_run_ac_pf"
        self.gridmodel.tell_solver_need_reset()
        expected_diff = 1e-2  # change slack impact is low
        load_p_init = 1.0 * np.array([el.target_p_mw for el in self.gridmodel.get_loads()])
        load_p_div = 5. * load_p_init
        
        # check that this load makes it diverge too
        tmp_grid = self.gridmodel.copy()
        [tmp_grid.change_p_load(l_id, val) for l_id, val in enumerate(load_p_div)]
        V_tmp =  tmp_grid.ac_pf(self.v_init, self.iter, self.tol_solver)
        assert len(V_tmp) == 0, "should have diverged !"
        
        self.aux_do_undo_ac(funname_do="_change_for_divergence_sbus",
                            funname_undo="_unchange_for_divergence_sbus",
                            runpf_fun=runpf_fun,
                            el_id=load_p_init,  # initial loads
                            expected_diff=expected_diff,
                            el_val=load_p_div,  # final loads
                            to_add_remove=0.,
                            )

    def _change_for_divergence_ybus(self, gridmodel, el_id, el_val):
        [gridmodel.deactivate_trafo(tr_id) for tr_id in el_id]
        
    def _unchange_for_divergence_ybus(self, gridmodel, el_id, el_val):
        [gridmodel.reactivate_trafo(tr_id) for tr_id in el_id]
            
    def test_divergence_ybus_ac(self, runpf_fun="_run_ac_pf"):
        """test I can make the grid diverge and converge again using ybus (islanding) in AC"""
        self.gridmodel.tell_solver_need_reset()
        expected_diff = 1e-2  # change slack impact is low
        trafo_to_disc = [0, 1, 2]
        
        # check that this load makes it diverge too
        tmp_grid = self.gridmodel.copy()
        [tmp_grid.deactivate_trafo(tr_id) for tr_id in trafo_to_disc]
        V_tmp =  tmp_grid.ac_pf(self.v_init, self.iter, self.tol_solver)
        assert len(V_tmp) == 0, "should have diverged !"
        
        self.aux_do_undo_ac(funname_do="_change_for_divergence_ybus",
                            funname_undo="_unchange_for_divergence_ybus",
                            runpf_fun=runpf_fun,
                            el_id=trafo_to_disc,  # trafo to disco
                            expected_diff=expected_diff,
                            el_val=0.,
                            to_add_remove=0.,
                            )

    def test_divergence_ybus_dc(self):
        """test I can make the grid diverge and converge again using ybus (islanding) in DC"""  
        self.test_divergence_ybus_ac("_run_dc_pf")   
        
    def _aux_disco_load(self, gridmodel, el_id, el_val):
        gridmodel.deactivate_load(el_id)
        
    def _aux_reco_load(self, gridmodel, el_id, el_val):
        gridmodel.reactivate_load(el_id)
        
    def test_disco_reco_load_ac(self, runpf_fun="_run_ac_pf"):
        """test I can disconnect a load (AC)"""
        for load in self.gridmodel.get_loads():
            self.gridmodel.tell_solver_need_reset()
            expected_diff = 1e-2
            if load.id == 3:
                expected_diff = 3e-3
            self.aux_do_undo_ac(funname_do="_aux_disco_load",
                                funname_undo="_aux_reco_load",
                                runpf_fun=runpf_fun,
                                el_id=load.id,
                                expected_diff=expected_diff,
                                el_val=0,
                                to_add_remove=0,
                                )

    def test_disco_reco_load_dc(self):   
        """test I can disconnect a load (DC)"""
        if not self.need_dc:
            self.skipTest("Useless to run DC")
        self.test_disco_reco_load_ac(runpf_fun="_run_dc_pf")
        
    def _aux_disco_gen(self, gridmodel, el_id, el_val):
        gridmodel.deactivate_gen(el_id)
        
    def _aux_reco_gen(self, gridmodel, el_id, el_val):
        gridmodel.reactivate_gen(el_id)
        
    def test_disco_reco_gen_ac(self, runpf_fun="_run_ac_pf"):
        """test I can disconnect a gen (AC)"""
        for gen in self.gridmodel.get_generators():
            if gen.is_slack:
                # by default single slack, so I don't disconnect it
                continue
            if gen.target_p_mw == 0.:
                # will have not impact as by default gen with p==0. are still pv
                continue
            self.gridmodel.tell_solver_need_reset()
            expected_diff = 1e-2
            if gen.id == 3:
                expected_diff = 3e-3
            self.aux_do_undo_ac(funname_do="_aux_disco_gen",
                                funname_undo="_aux_reco_gen",
                                runpf_fun=runpf_fun,
                                el_id=gen.id,
                                expected_diff=expected_diff,
                                el_val=0,
                                to_add_remove=0,
                                )

    def test_disco_reco_gen_dc(self):   
        """test I can disconnect a shunt (DC)"""
        if not self.need_dc:
            self.skipTest("Useless to run DC")
        self.test_disco_reco_gen_ac(runpf_fun="_run_dc_pf")
        
    def _aux_disco_shunt(self, gridmodel, el_id, el_val):
        gridmodel.deactivate_shunt(el_id)
        
    def _aux_reco_shunt(self, gridmodel, el_id, el_val):
        gridmodel.reactivate_shunt(el_id)
        
    def test_disco_reco_shunt_ac(self, runpf_fun="_run_ac_pf"):
        """test I can disconnect a shunt (AC)"""
        for shunt in self.gridmodel.get_shunts():
            self.gridmodel.tell_solver_need_reset()
            expected_diff = 1e-2
            if runpf_fun == "_run_dc_pf":
                # in dc q is not used, so i skipped if no target_p_mw
                if shunt.target_p_mw == 0:
                    continue
            self.aux_do_undo_ac(funname_do="_aux_disco_shunt",
                                funname_undo="_aux_reco_shunt",
                                runpf_fun=runpf_fun,
                                el_id=shunt.id,
                                expected_diff=expected_diff,
                                el_val=0,
                                to_add_remove=0,
                                )

    def test_disco_reco_shunt_dc(self):   
        """test I can disconnect a shunt (DC)"""
        if not self.need_dc:
            self.skipTest("Useless to run DC")
        self.test_disco_reco_shunt_ac(runpf_fun="_run_dc_pf")
    
    def _aux_disco_shunt(self, gridmodel, el_id, el_val):
        gridmodel.deactivate_bus(0)
        gridmodel.reactivate_bus(15)
        gridmodel.change_bus_powerline_or(0, 15)
        gridmodel.change_bus_powerline_or(1, 15)
        gridmodel.change_bus_gen(5, 15)
        
    def _aux_reco_shunt(self, gridmodel, el_id, el_val):
        gridmodel.reactivate_bus(0)
        gridmodel.deactivate_bus(15)
        gridmodel.change_bus_powerline_or(0, 0)
        gridmodel.change_bus_powerline_or(1, 0)
        gridmodel.change_bus_gen(5, 0)
        
    def test_change_bus2_ac(self, runpf_fun="_run_ac_pf"):
        """test for bus 2, basic test I don't do it for all kind of objects (AC pf)"""
        expected_diff = 1e-2
        self.aux_do_undo_ac(funname_do="_aux_disco_shunt",
                            funname_undo="_aux_reco_shunt",
                            runpf_fun=runpf_fun,
                            el_id=0,
                            expected_diff=expected_diff,
                            el_val=0,
                            to_add_remove=0,
                            )
        
 
class TestSolverControlNRSing(TestSolverControl):
    def _aux_setup_grid(self):
        self.need_dc = False  # is it worth it to run DC powerflow ?
        self.can_dist_slack = False
        self.gridmodel.change_solver(SolverType.SparseLUSingleSlack)
        self.gridmodel.change_solver(SolverType.DC)
 
 
class TestSolverControlFDPF_XB(TestSolverControl):
    def _aux_setup_grid(self):
        self.need_dc = False  # is it worth it to run DC powerflow ?
        self.can_dist_slack = False
        self.gridmodel.change_solver(SolverType.FDPF_XB_SparseLU)
        self.gridmodel.change_solver(SolverType.DC)
        self.iter = 30
 
 
class TestSolverControlFDPF_BX(TestSolverControl):
    def _aux_setup_grid(self):
        self.need_dc = False  # is it worth it to run DC powerflow ?
        self.can_dist_slack = False
        self.gridmodel.change_solver(SolverType.FDPF_BX_SparseLU)
        self.gridmodel.change_solver(SolverType.DC)
        self.iter = 30
               
 
class TestSolverControlGaussSeidel(TestSolverControl):
    def _aux_setup_grid(self):
        self.need_dc = False  # is it worth it to run DC powerflow ?
        self.can_dist_slack = False
        self.gridmodel.change_solver(SolverType.GaussSeidel)
        self.iter = 350
 
 
class TestSolverControlGaussSeidelSynch(TestSolverControl):
    def _aux_setup_grid(self):
        self.need_dc = False  # is it worth it to run DC powerflow ?
        self.can_dist_slack = False
        self.gridmodel.change_solver(SolverType.GaussSeidelSynch)
        self.iter = 1000
               
               
if __name__ == "__main__":
    unittest.main()
