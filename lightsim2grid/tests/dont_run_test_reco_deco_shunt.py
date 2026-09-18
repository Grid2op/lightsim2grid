import sys
import grid2op
from grid2op.tests.BaseBackendTest import (BaseTestNames,  # noqa: E402
                                           BaseTestLoadingCase,
                                           BaseTestLoadingBackendFunc, 
                                           BaseTestTopoAction,
                                           BaseTestEnvPerformsCorrectCascadingFailures,
                                           BaseTestChangeBusAffectRightBus,
                                        #    BaseTestShuntAction,
                                           BaseTestResetEqualsLoadGrid,
                                           BaseTestVoltageOWhenDisco,
                                           BaseTestChangeBusSlack,
                                           BaseIssuesTest,
                                           BaseStatusActions,
                                           MakeBackend, AlwaysLegal, CompleteAction, AmbiguousAction
                                           )

# gdb python --eval-command "run test_reco_deco_shunt.py"
from lightsim2grid import LightSimBackend

backend2 = LightSimBackend()
print("MAKE HERE")
env_change_q = grid2op.make(
    "rte_case14_realistic",
    test=True,
    gamerules_class=AlwaysLegal,
    action_class=CompleteAction,
    backend=backend2,
    _add_to_name="_3"
)
assert not env_change_q.done
print("MAKE END")

param = env_change_q.parameters
param.NO_OVERFLOW_DISCONNECTION = True
env_change_q.change_parameters(param)
print("RESET HERE")
env_change_q.reset(seed=0, options={"time serie id": 0})
print("RESET END")

# obs_disco_sh, *_ = env_change_q.step(
#             env_change_q.action_space({"shunt": {"set_bus": [(0, -1)]}})
#         )
# obs_disco_sh, *_ = env_change_q.step(
#             env_change_q.action_space({"shunt": {"set_bus": [(0, -1)]}})
#         )


obs_co_bus2_sh_alone, reward, done, info = env_change_q.step(
    env_change_q.action_space({"shunt": {"set_bus": [(0, 2)]}})
)
print(f"{done = }")
env_change_q.reset(seed=0, options={"time serie id": 0})
