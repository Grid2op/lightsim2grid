# this class is defined to avoid circular references when using grid2op -> pandapower -> lightsim2grid -> grid2op
# this module is lazily imported by the LightSimBackend and should not be used anywhere else, 
# as the name states

import warnings

try:
    import grid2op
    from grid2op.Backend import PandaPowerBackend
    class _DoNotUseAnywherePandaPowerBackend(PandaPowerBackend):
        """used to duplicate the class attributes of PandaPowerBackend"""
        shunts_data_available = True
        
    if not hasattr(PandaPowerBackend, "_clear_grid_dependant_class_attributes"):
        # not available in python 3.7
        warnings.warn("This is a compatibility mode, please avoid using it and migrate to python >=3.10 ideally.")
        def f(cls):
            cls._INIT_GRID_CLS = None  # do not modify that, this is handled by grid2op automatically
            cls._PATH_GRID_CLASSES = None  # especially do not modify that
            
            cls.glop_version = grid2op.__version__

            cls.SUB_COL = 0
            cls.LOA_COL = 1
            cls.GEN_COL = 2
            cls.LOR_COL = 3
            cls.LEX_COL = 4
            cls.STORAGE_COL = 5

            cls.attr_list_vect = None
            cls.attr_list_set = {}
            cls.attr_list_json = []
            cls.attr_nan_list_set = set()

            # class been init
            cls._IS_INIT = False

            # name of the objects
            cls.env_name = "unknown"
            cls.name_load = None
            cls.name_gen = None
            cls.name_line = None
            cls.name_sub = None
            cls.name_storage = None

            cls.n_gen = -1
            cls.n_load = -1
            cls.n_line = -1
            cls.n_sub = -1
            cls.n_storage = -1

            cls.sub_info = None
            cls.dim_topo = -1

            # to which substation is connected each element
            cls.load_to_subid = None
            cls.gen_to_subid = None
            cls.line_or_to_subid = None
            cls.line_ex_to_subid = None
            cls.storage_to_subid = None

            # which index has this element in the substation vector
            cls.load_to_sub_pos = None
            cls.gen_to_sub_pos = None
            cls.line_or_to_sub_pos = None
            cls.line_ex_to_sub_pos = None
            cls.storage_to_sub_pos = None

            # which index has this element in the topology vector
            cls.load_pos_topo_vect = None
            cls.gen_pos_topo_vect = None
            cls.line_or_pos_topo_vect = None
            cls.line_ex_pos_topo_vect = None
            cls.storage_pos_topo_vect = None

            # "convenient" way to retrieve information of the grid
            cls.grid_objects_types = None
            # to which substation each element of the topovect is connected
            cls._topo_vect_to_sub = None

            # list of attribute to convert it from/to a vector
            cls._vectorized = None

            # redispatch data, not available in all environment
            cls.redispatching_unit_commitment_availble = False
            cls.gen_type = None
            cls.gen_pmin = None
            cls.gen_pmax = None
            cls.gen_redispatchable = None
            cls.gen_max_ramp_up = None
            cls.gen_max_ramp_down = None
            cls.gen_min_uptime = None
            cls.gen_min_downtime = None
            cls.gen_cost_per_MW = None  # marginal cost (in currency / (power.step) and not in $/(MW.h) it would be $ / (MW.5mins) )
            cls.gen_startup_cost = None  # start cost (in currency)
            cls.gen_shutdown_cost = None  # shutdown cost (in currency)
            cls.gen_renewable = None

            # storage unit static data
            cls.storage_type = None
            cls.storage_Emax = None
            cls.storage_Emin = None
            cls.storage_max_p_prod = None
            cls.storage_max_p_absorb = None
            cls.storage_marginal_cost = None
            cls.storage_loss = None
            cls.storage_charging_efficiency = None
            cls.storage_discharging_efficiency = None

            # grid layout
            cls.grid_layout = None

            # shunt data, not available in every backend
            cls.n_shunt = None
            cls.name_shunt = None
            cls.shunt_to_subid = None

            # alarm / alert
            cls.assistant_warning_type = None
            
            # alarms
            cls.dim_alarms = 0
            cls.alarms_area_names = []
            cls.alarms_lines_area = {}
            cls.alarms_area_lines = []

            # alerts
            cls.dim_alerts = 0
            cls.alertable_line_names = []
            cls.alertable_line_ids = []
        setattr(_DoNotUseAnywherePandaPowerBackend, "_clear_grid_dependant_class_attributes", classmethod(f))
            
        
except ImportError as exc_:
    # grid2op is not installed, we do not use it.
    pass
