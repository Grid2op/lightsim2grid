# this class is defined to avoid circular references when using grid2op -> pandapower -> lightsim2grid -> grid2op
# this module is lazily imported by the LightSimBackend and should not be used anywhere else, 
# as the name states
try:
    from grid2op.Backend import PandaPowerBackend
    class _DoNotUseAnywherePandaPowerBackend(PandaPowerBackend):
        """used to duplicate the class attributes of PandaPowerBackend"""
        pass
except ImportError as exc_:
    # grid2op is not installed, we do not use it.
    pass