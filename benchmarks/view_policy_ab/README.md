# Result views: `reference` vs `reference_internal`

The A/B behind `docs/dev_notes/result_view_lifetime.md`: what binding the result getters of
`LSGrid` (`get_gen_res`, `get_line_res1_full`, ...) with `py::return_value_policy::reference_internal`
instead of `reference` costs per call.

`ab_views.cpp` binds the very same C++ getters three times -- `reference` (what
`binding_lsgrid.cpp` did before 1.3.0), `reference_internal` (each array's `.base` is the grid,
which it then keeps alive: what it does since) and `copy` (each array owns a copy, for scale) --
on a synthetic meshed grid of any size. Nothing else differs: same pybind11, same casters, same
object, same process.

```
python -m pip install pybind11 numpy
cmake -S benchmarks/view_policy_ab -B build_ab -DCMAKE_BUILD_TYPE=Release \
      -Dpybind11_DIR=$(python -m pybind11 --cmakedir)
cmake --build build_ab -j
AB_VIEWS_DIR=build_ab python benchmarks/view_policy_ab/run_ab.py time   # ns per call
AB_VIEWS_DIR=build_ab python benchmarks/view_policy_ab/run_ab.py ir     # instructions per call (valgrind)
```

`time` interleaves the variants in shuffled rounds and reports the median and the minimum per
call, with a no-op method for the dispatch alone. `ir` runs each variant under callgrind twice
with a different number of calls and takes the difference, so that the interpreter start, the
imports and the grid build cancel out; the counts are deterministic.

Note what a per-call cost means for grid2op: `LightSimBackend._fetch_grid_data` fetches each
result view once and reads it again after every powerflow, so the backend pays it once per view,
not per step.
