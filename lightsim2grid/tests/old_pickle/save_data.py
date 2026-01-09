import os
import pickle
import warnings
import grid2op
from lightsim2grid import LightSimBackend


def main():    
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        env = grid2op.make("l2rpn_idf_2023", test=True, backend=LightSimBackend())    
    gridmodel = env.backend._grid
    
    with open(os.path.join("test_pickle.pickle"), "wb") as f:
        pickle.dump(gridmodel, f)
        
    for fun_name in [
        "get_loads",
        "get_lines",
        "get_trafos",
        "get_storages",
        "get_generators",
        "get_shunts",
        "get_substations",
        ]:
        with open(os.path.join(f"test_pickle_{fun_name}.pickle"), "wb") as f:
                pickle.dump(getattr(gridmodel, fun_name)(), f)
                
if __name__ == "__main__":
    main()
