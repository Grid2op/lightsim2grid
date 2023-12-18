# TODO: test that "if I do something, then with or without the speed optim, I have the same thing"

# use backend._grid.tell_solver_need_reset() to force the reset of the solver
# and by "something" do :
# - disconnect line X
# - reconnect line X
# - change load
# - change gen
# - change shunt
# - change slack
# - change slack weight
# - when it diverges (and then I can make it un converge normally)
# - test change bus to -1 and deactivate element has the same impact

# TODO and do that for all solver type: NR, NRSingleSlack, GaussSeidel, GaussSeidelSynch, FDPF_XB and FDPFBX
# and for all solver check everything that can be check: final tolerance, number of iteration, etc. etc.

import unittest
import warnings
import numpy as np
import grid2op
from grid2op.Action import CompleteAction
from lightsim2grid import LightSimBackend


class TestSolverControl(unittest.TestCase):
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("educ_case14_storage",
                                    test=True,
                                    action_class=CompleteAction,
                                    backend=LightSimBackend())
        self.gridmodel = self.env.backend._grid
        self.v_init = 1.0 * self.env.backend.V
        self.iter = 10
        self.tol_solver = 1e-8  # solver
        self.tol_equal = 1e-10  #  for comparing with and without the "smarter solver" things, and make sure everything is really equal!
    
    def test_pf_run_dc(self):
        """test I have the same results if nothing is done with and without restarting from scratch when running dc powerflow"""
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
                       expected_diff=0.1
                       ):
        pf_mode = "AC" if runpf_fun=="_run_ac_pf" else "DC"
        print(f"Looking at {el_id} in {pf_mode}")
        
        # test "do the action"
        V_init = getattr(self, runpf_fun)(gridmodel=self.gridmodel)
        assert len(V_init), f"error for el_id={el_id}: gridmodel should converge in {pf_mode}"
        
        getattr(self, funname_do)(gridmodel=self.gridmodel, el_id=el_id, el_val=el_val)
        V_disc = getattr(self, runpf_fun)(gridmodel=self.gridmodel)
        if len(V_disc) > 0:
            # powerflow converges, all should converge
            assert (np.abs(V_init - V_disc) >= expected_diff).any(), f"error for el_id={el_id}: at least one bus should have changed its result voltage in {pf_mode}"
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
        assert np.allclose(V_reco, V_init, rtol=self.tol_equal, atol=self.tol_equal), f"error for el_id={el_id}: do an action and then undo it should not have any impact in {pf_mode}"
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
        self.test_disco_reco_trafo_ac(runpf_fun="_run_dc_pf")
            