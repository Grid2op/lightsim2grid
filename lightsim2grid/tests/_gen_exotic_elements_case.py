# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Dev-only generator: captures the `LSGrid.init_*`/`add_*`/`set_*` calls made by
`init_from_pypowsybl` while building `_exotic_elements_fixture.build_exotic_elements_grid()`,
and emits two "case files" that reproduce the same grid via hardcoded literals, with no
pypowsybl/pandapower/grid2op dependency at runtime:

* `lightsim2grid/tests/_exotic_elements_case_ls.py` (numpy + lightsim2grid only)
* `src/tests/case_exotic_elements.hpp` (pure C++, Eigen + lightsim2grid core only)

How it works: `LSGrid` (`lightsim2grid.network.LSGrid`) is a plain, mutable pybind11
class, so its grid-construction methods can be monkeypatched at the class level for the
duration of a single real `init_from_pypowsybl(...)` call -- every allowlisted method
records its (name, args) before delegating to the real implementation, and the patch is
reverted immediately after. This means the actual `initLSGrid.init()` orchestration is
never duplicated here: if it changes, this generator picks up the new call sequence
automatically the next time it is rerun.

Only calls that affect the solve are kept (the allowlist below); purely cosmetic /
bookkeeping calls (`set_*_names`, `set_bus_voltage_limits`, `set_*_current_limit_*`,
`set_*_to_subid`) are recognized and dropped -- matching how `src/tests/test_lsgrid.cpp`'s
own hand-written helpers never touch any of those either. Any call observed that is
neither allowlisted nor in the explicit skip set raises immediately, so a future change to
`initLSGrid.py` cannot silently produce an incomplete case file.

Regenerate with (requires pypowsybl; NOT part of CI, NOT imported by any test module):

    venv_ls/bin/python -m lightsim2grid.tests._gen_exotic_elements_case

...whenever `_exotic_elements_fixture.build_exotic_elements_network()` changes.
"""

import datetime
import os

import numpy as np
import pypowsybl

import lightsim2grid
from lightsim2grid.network import LSGrid, init_from_pypowsybl

from ._exotic_elements_fixture import build_exotic_elements_network

# ---------------------------------------------------------------------------
# 1. call capture
# ---------------------------------------------------------------------------

# Calls that affect the AC/DC solve: kept, emitted verbatim into both case files.
ALLOWLIST = {
    "set_sn_mva", "set_init_vm_pu",
    "init_bus",
    "init_generators_full", "deactivate_gen", "reactivate_gen", "set_gen_regulated_bus",
    "init_loads", "deactivate_load", "reactivate_load",
    "init_powerlines_full",
    "deactivate_powerline", "deactivate_powerline_side1", "deactivate_powerline_side2",
    "reactivate_powerline",
    "init_trafo",
    "deactivate_trafo", "deactivate_trafo_side1", "deactivate_trafo_side2",
    "reactivate_trafo",
    "set_trafo_shift_dependent_rx",
    "init_shunt", "deactivate_shunt", "reactivate_shunt",
    "init_svcs", "deactivate_svc", "reactivate_svc",
    "init_hvdc_lines",
    "deactivate_dcline", "deactivate_dcline_side1", "deactivate_dcline_side2",
    "reactivate_dcline",
    "init_storages", "deactivate_storage", "reactivate_storage",
    "add_gen_slackbus", "remove_gen_slackbus",
    "consider_only_main_component", "init_bus_status",
}

# Purely cosmetic / grid2op-bookkeeping calls: recognized and silently dropped.
SKIP_SET = {
    "set_substation_names", "set_bus_voltage_limits",
    "set_gen_names", "set_load_names",
    "set_line_names", "set_line_current_limit_side1", "set_line_current_limit_side2",
    "set_trafo_names", "set_trafo_current_limit_side1", "set_trafo_current_limit_side2",
    "set_shunt_names", "set_svc_names", "set_dcline_names", "set_storage_names",
    "set_gen_to_subid", "set_load_to_subid", "set_storage_to_subid", "set_shunt_to_subid",
    "set_line_to_sub1_id", "set_line_to_sub2_id", "set_trafo_to_sub1_id", "set_trafo_to_sub2_id",
    "tell_solver_need_reset",
}


class _RecordingLSGrid:
    """Context manager: monkeypatches every `ALLOWLIST`/`SKIP_SET` method found on the
    `LSGrid` class to log `(name, args)` (kwargs are asserted empty -- the C++ bindings
    are positional-only) before delegating to the real implementation, and restores the
    originals on exit. `self.call_log` only ever contains allowlisted calls."""

    def __init__(self):
        self.call_log = []
        self._originals = {}

    def __enter__(self):
        for name in ALLOWLIST | SKIP_SET:
            if hasattr(LSGrid, name):
                self._originals[name] = getattr(LSGrid, name)
        for name, original in self._originals.items():
            setattr(LSGrid, name, self._make_wrapper(name, original))
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        for name, original in self._originals.items():
            setattr(LSGrid, name, original)
        return False

    def _make_wrapper(self, name, original):
        keep = name in ALLOWLIST

        def wrapper(this, *args, **kwargs):
            assert not kwargs, (
                f"{name}: unexpected kwargs {kwargs!r} -- the generator assumes every "
                "LSGrid call from init_from_pypowsybl is positional-only")
            if keep:
                self.call_log.append((name, args))
            return original(this, *args)

        return wrapper


def _capture_exotic_elements_calls():
    net = build_exotic_elements_network()
    with _RecordingLSGrid() as recorder:
        init_from_pypowsybl(net, gen_slack_id="B1-G")
    seen_names = {name for name, _ in recorder.call_log}
    unexpected = seen_names - ALLOWLIST
    assert not unexpected, f"unclassified LSGrid calls observed: {unexpected}"
    return recorder.call_log


# ---------------------------------------------------------------------------
# 2. value formatting shared by both languages
# ---------------------------------------------------------------------------

def _double_literal(x):
    """`repr()` on a python float is the shortest decimal string that round-trips to
    the exact same IEEE-754 double, and that string is valid literal syntax in both
    Python and C++ -- so a single helper backs both formatters and guarantees the
    regenerated grid reaches the bit-identical Newton-Raphson fixed point."""
    x = float(x)
    assert np.isfinite(x), f"non-finite value {x!r} encountered; formatter does not handle nan/inf"
    return repr(x)


def _value_kind(value):
    """Classify a captured argument into one of the shapes the LSGrid C++ API accepts,
    driving both the Python and the C++ formatter."""
    if isinstance(value, (bool, np.bool_)):
        return "scalar_bool"
    if isinstance(value, (int, np.integer)):
        return "scalar_int"
    if isinstance(value, (float, np.floating)):
        return "scalar_real"
    if isinstance(value, np.ndarray):
        if np.issubdtype(value.dtype, np.complexfloating):
            return "cplx_vect"
        if np.issubdtype(value.dtype, np.bool_):
            return "vector_bool"
        if np.issubdtype(value.dtype, np.integer):
            return "int_vect"
        if np.issubdtype(value.dtype, np.floating):
            return "real_vect"
        raise TypeError(f"unhandled ndarray dtype {value.dtype}")
    if isinstance(value, list):
        assert len(value) > 0, "empty python list argument: ambiguous element type, not handled"
        if isinstance(value[0], list):
            return "nested_real"
        if isinstance(value[0], bool):
            return "vector_bool"
        if isinstance(value[0], int):
            return "vector_int"
        raise TypeError(f"unhandled list element type {type(value[0])}")
    raise TypeError(f"unhandled argument type {type(value)}")


def _py_format_value(value):
    kind = _value_kind(value)
    if kind == "scalar_bool":
        return repr(bool(value))
    if kind == "scalar_int":
        return repr(int(value))
    if kind == "scalar_real":
        return _double_literal(value)
    if kind == "real_vect":
        elems = ", ".join(_double_literal(v) for v in value)
        return f"np.array([{elems}], dtype=np.{value.dtype.name})"
    if kind == "int_vect":
        elems = ", ".join(str(int(v)) for v in value)
        return f"np.array([{elems}], dtype=np.{value.dtype.name})"
    if kind == "cplx_vect":
        elems = ", ".join(f"complex({_double_literal(v.real)}, {_double_literal(v.imag)})" for v in value)
        return f"np.array([{elems}], dtype=np.{value.dtype.name})"
    if kind == "vector_bool":
        return "[" + ", ".join(repr(bool(v)) for v in value) + "]"
    if kind == "vector_int":
        return "[" + ", ".join(repr(int(v)) for v in value) + "]"
    if kind == "nested_real":
        rows = ("[" + ", ".join(_double_literal(v) for v in row) + "]" for row in value)
        return "[" + ", ".join(rows) + "]"
    raise AssertionError(kind)  # pragma: no cover


def _cpp_format_arg(varname, value):
    """Returns (decl_lines, expr): zero or more C++ statements to declare/fill `varname`,
    and the expression to use for this argument in the call itself."""
    kind = _value_kind(value)
    if kind == "scalar_bool":
        return [], ("true" if value else "false")
    if kind == "scalar_int":
        return [], str(int(value))
    if kind == "scalar_real":
        return [], _double_literal(value)
    if kind == "real_vect":
        n = len(value)
        decl = [f"RealVect {varname}({n});"]
        if n > 0:
            decl.append(f"{varname} << " + ", ".join(_double_literal(v) for v in value) + ";")
        return decl, varname
    if kind == "int_vect":
        n = len(value)
        decl = [f"Eigen::VectorXi {varname}({n});"]
        if n > 0:
            decl.append(f"{varname} << " + ", ".join(str(int(v)) for v in value) + ";")
        return decl, varname
    if kind == "cplx_vect":
        n = len(value)
        decl = [f"CplxVect {varname}({n});"]
        if n > 0:
            elems = ", ".join(
                f"std::complex<double>({_double_literal(v.real)}, {_double_literal(v.imag)})" for v in value)
            decl.append(f"{varname} << {elems};")
        return decl, varname
    if kind == "vector_bool":
        elems = ", ".join("true" if v else "false" for v in value)
        return [f"const std::vector<bool> {varname}{{{elems}}};"], varname
    if kind == "vector_int":
        elems = ", ".join(str(int(v)) for v in value)
        return [f"const std::vector<int> {varname}{{{elems}}};"], varname
    if kind == "nested_real":
        decl = [f"const std::vector<std::vector<real_type>> {varname} = {{"]
        for row in value:
            row_lits = ", ".join(_double_literal(v) for v in row)
            decl.append(f"    {{{row_lits}}},")
        decl.append("};")
        return decl, varname
    raise AssertionError(kind)  # pragma: no cover


# Readable C++ variable names, taken straight from the LSGrid.hpp parameter names, for
# every allowlisted method with array-typed arguments. Methods not listed here (all
# scalar args, e.g. the deactivate_* family) fall back to `f"{method}_arg{i}"`.
CPP_PARAM_NAMES = {
    "init_bus": ["n_sub", "n_busbar_per_sub", "bus_vn_kv", "nb_line", "nb_trafo"],
    "init_generators_full": [
        "gen_p", "gen_v_pu", "gen_q", "gen_voltage_regulator_on", "gen_min_q", "gen_max_q", "gen_bus_id"],
    "init_loads": ["load_p", "load_q", "load_bus_id"],
    "init_powerlines_full": ["line_r", "line_x", "line_h1", "line_h2", "line_from_id", "line_to_id"],
    "init_trafo": [
        "trafo_r", "trafo_x", "trafo_b", "trafo_ratio", "trafo_shift_degree",
        "trafo_tap_hv", "trafo_bus1_id", "trafo_bus2_id", "trafo_ignore_tap_side_for_shift"],
    "set_trafo_shift_dependent_rx": ["trafo_shift_rx_enable", "trafo_shift_alpha_rad", "trafo_shift_rx_corr_pct"],
    "init_shunt": ["shunt_p_mw", "shunt_q_mvar", "shunt_bus_id"],
    "init_svcs": [
        "svc_regulation_mode", "svc_target_vm_pu", "svc_q_setpoint_mvar", "svc_slope_pu",
        "svc_b_min", "svc_b_max", "svc_regulated_bus_id", "svc_bus_id"],
    "init_hvdc_lines": [
        "hvdc_bus1_id", "hvdc_bus2_id", "hvdc_type1", "hvdc_type2",
        "hvdc_loss_factor1", "hvdc_loss_factor2", "hvdc_vreg1_on", "hvdc_vreg2_on",
        "hvdc_vm1_pu", "hvdc_vm2_pu", "hvdc_q1_setpoint_mvar", "hvdc_q2_setpoint_mvar",
        "hvdc_min_q1", "hvdc_max_q1", "hvdc_min_q2", "hvdc_max_q2",
        "hvdc_power_factor1", "hvdc_power_factor2", "hvdc_converters_mode",
        "hvdc_p_setpoint_mw", "hvdc_r_ohm", "hvdc_nominal_v_kv",
        "hvdc_droop_enabled", "hvdc_droop_p0_mw", "hvdc_droop_mw_per_deg",
        "hvdc_pmax_1to2_mw", "hvdc_pmax_2to1_mw"],
    "init_storages": ["storage_p", "storage_q", "storage_bus_id"],
}


# ---------------------------------------------------------------------------
# 3. rendering
# ---------------------------------------------------------------------------

_LICENSE_PY = '''\
# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.
'''

_LICENSE_CPP = '''\
// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.
'''

_REGEN_COMMAND = "python -m lightsim2grid.tests._gen_exotic_elements_case"


def _generation_stamp():
    """One line, reused verbatim (same day, same tool versions) in both generated
    files, so a diff of either one always shows when, how, and against which
    lightsim2grid/pypowsybl versions it was last regenerated."""
    today = datetime.date.today().isoformat()
    ls_version = lightsim2grid.__version__
    pp_version = pypowsybl.__version__
    return (f"Generated on {today} using `{_REGEN_COMMAND}` "
            f"(lightsim2grid {ls_version}, pypowsybl {pp_version}).")


def render_python(call_log, stamp):
    lines = [_LICENSE_PY]
    lines.append('"""AUTO-GENERATED by lightsim2grid/tests/_gen_exotic_elements_case.py -- DO NOT EDIT BY HAND.')
    lines.append('')
    lines.append('Reproduces the grid built by `_exotic_elements_fixture.build_exotic_elements_grid()` (IEEE14 +')
    lines.append('SVC + storage + 3 HVDC lines + a phase-shifting transformer) via literal `LSGrid.init_*`/')
    lines.append('`add_*`/`set_*` calls captured from a real `init_from_pypowsybl` run -- no pypowsybl/')
    lines.append('pandapower/grid2op dependency at runtime, only numpy + lightsim2grid.')
    lines.append('')
    lines.append(stamp)
    lines.append('"""')
    lines.append('')
    lines.append('import numpy as np')
    lines.append('')
    lines.append('from lightsim2grid.network import LSGrid')
    lines.append('')
    lines.append('')
    lines.append('def build_exotic_elements_case_grid() -> LSGrid:')
    lines.append('    """Unsolved: callers run their own ac_pf / dc_pf, same contract as')
    lines.append('    `_exotic_elements_fixture.build_exotic_elements_grid(solve=False)`."""')
    lines.append('    model = LSGrid()')
    lines.append('')
    for name, args in call_log:
        arg_strs = [_py_format_value(a) for a in args]
        if not arg_strs:
            lines.append(f"    model.{name}()")
            continue
        lines.append(f"    model.{name}(")
        for i, a in enumerate(arg_strs):
            comma = "," if i < len(arg_strs) - 1 else ""
            lines.append(f"        {a}{comma}")
        lines.append("    )")
    lines.append('')
    lines.append('    return model')
    lines.append('')
    return "\n".join(lines)


def render_cpp(call_log, stamp):
    lines = [_LICENSE_CPP]
    lines.append('#ifndef LS2G_CASE_EXOTIC_ELEMENTS_H')
    lines.append('#define LS2G_CASE_EXOTIC_ELEMENTS_H')
    lines.append('')
    lines.append('// AUTO-GENERATED by lightsim2grid/tests/_gen_exotic_elements_case.py -- DO NOT EDIT BY HAND.')
    lines.append('//')
    lines.append('// Reproduces the grid built by lightsim2grid.tests._exotic_elements_fixture.')
    lines.append('// build_exotic_elements_grid() (IEEE14 + SVC + storage + 3 HVDC lines + a phase-shifting')
    lines.append('// transformer) via literal LSGrid init_*/add_*/set_* calls captured from a real')
    lines.append('// init_from_pypowsybl run -- no pypowsybl/pandapower/grid2op dependency, pure C++.')
    lines.append('//')
    lines.append(f'// {stamp}')
    lines.append('')
    lines.append('#include <complex>')
    lines.append('#include <vector>')
    lines.append('')
    lines.append('#include "LSGrid.hpp"')
    lines.append('')
    lines.append('namespace ls2g_test {')
    lines.append('')
    lines.append('inline ls2g::LSGrid make_exotic_elements_grid()')
    lines.append('{')
    lines.append('    using ls2g::LSGrid;')
    lines.append('    using ls2g::RealVect;')
    lines.append('    using ls2g::CplxVect;')
    lines.append('    using ls2g::real_type;')
    lines.append('')
    lines.append('    LSGrid grid;')
    lines.append('')
    for name, args in call_log:
        param_names = CPP_PARAM_NAMES.get(name)
        call_arg_exprs = []
        for i, a in enumerate(args):
            varname = param_names[i] if param_names and i < len(param_names) else f"{name}_arg{i}"
            decl_lines, expr = _cpp_format_arg(varname, a)
            for decl_line in decl_lines:
                lines.append(f"    {decl_line}")
            call_arg_exprs.append(expr)
        if call_arg_exprs:
            lines.append(f"    grid.{name}(" + ", ".join(call_arg_exprs) + ");")
        else:
            lines.append(f"    grid.{name}();")
        lines.append('')
    lines.append('    return grid;')
    lines.append('}')
    lines.append('')
    lines.append('}  // namespace ls2g_test')
    lines.append('')
    lines.append('#endif  // LS2G_CASE_EXOTIC_ELEMENTS_H')
    lines.append('')
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 4. entry point
# ---------------------------------------------------------------------------

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PY_CASE_PATH = os.path.join(_REPO_ROOT, "lightsim2grid", "tests", "_exotic_elements_case_ls.py")
CPP_CASE_PATH = os.path.join(_REPO_ROOT, "src", "tests", "case_exotic_elements.hpp")


def main():
    call_log = _capture_exotic_elements_calls()
    stamp = _generation_stamp()

    py_src = render_python(call_log, stamp)
    with open(PY_CASE_PATH, "w") as f:
        f.write(py_src)
    print(f"wrote {PY_CASE_PATH} ({len(call_log)} calls)")

    cpp_src = render_cpp(call_log, stamp)
    with open(CPP_CASE_PATH, "w") as f:
        f.write(cpp_src)
    print(f"wrote {CPP_CASE_PATH} ({len(call_log)} calls)")


if __name__ == "__main__":
    main()
