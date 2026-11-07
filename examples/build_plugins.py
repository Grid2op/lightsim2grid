#!/usr/bin/env python3
"""Configure and build the example solver plugins in one go, on Linux, macOS and Windows.

Every sub-folder of ``examples/`` holding a ``CMakeLists.txt`` is a plugin. Each one is built
in its own ``<plugin>/build`` folder, which is where its ``test_plugin.py`` looks for it, and
where ``lightsim2grid.load_algorithm_plugin("<plugin>")`` finds it by name when lightsim2grid
is an editable install of this checkout.

By default the plugins are built against the *installed* lightsim2grid (the one
``python -c "import lightsim2grid"`` finds), whose CMake package directory is asked from a
scratch directory so that the source tree is not imported instead.

Command line::

    python examples/build_plugins.py                        # every plugin
    python examples/build_plugins.py rebuild_outer lm       # a subset ("_algorithm" optional)
    python examples/build_plugins.py --list
    python examples/build_plugins.py -j 4 --test            # build, then run each test_plugin.py
    python examples/build_plugins.py --no-clean             # incremental rebuild
    python examples/build_plugins.py --build-type Debug -D CMAKE_VERBOSE_MAKEFILE=ON
    python examples/build_plugins.py --python .venv/bin/python
    python examples/build_plugins.py --source-tree          # headers of this checkout, no install

From Python::

    from build_plugins import build_plugins
    build_plugins(["rebuild_outer_algorithm"], jobs=4, test=True)
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

EXAMPLES_DIR = Path(__file__).resolve().parent
_NOT_PLUGINS = {"cmake"}  # examples/cmake holds a helper included by the plugins
_LIB_SUFFIXES = (".so", ".dll", ".dylib")


def available_plugins():
    """The plugin folders under examples/, sorted by name."""
    return sorted(p.name for p in EXAMPLES_DIR.iterdir()
                  if p.is_dir() and p.name not in _NOT_PLUGINS and (p / "CMakeLists.txt").is_file())


def _resolve_names(names):
    """Accept the folder name, or the name without its "_algorithm" suffix."""
    all_plugins = available_plugins()
    if not names:
        return all_plugins
    res = []
    for name in names:
        name = name.rstrip("/\\").split("/")[-1].split("\\")[-1]
        for cand in (name, f"{name}_algorithm"):
            if cand in all_plugins:
                res.append(cand)
                break
        else:
            raise ValueError(f"Unknown plugin '{name}'. Available: {', '.join(all_plugins)}")
    return res


def get_installed_cmake_dir(python=None):
    """lightsim2grid.get_cmake_dir() of the installed package.

    Asked from a temporary directory: run from the repository root, ``import lightsim2grid``
    would pick the source tree, which has no CMake package. When the import fails (an extension
    linked to a shared library not on the loader path, e.g. NICSLU / CKTSO), the package is
    looked for in the interpreter's site-packages instead, which loads nothing.
    """
    python = python or sys.executable
    with tempfile.TemporaryDirectory() as tmp:
        proc = subprocess.run(
            [python, "-c", "import lightsim2grid; print(lightsim2grid.get_cmake_dir())"],
            cwd=tmp, capture_output=True, text=True)
        if proc.returncode == 0:
            return Path(proc.stdout.strip().splitlines()[-1])
        fallback = subprocess.run(
            [python, "-c",
             "import os, sysconfig\n"
             "p = sysconfig.get_paths()\n"
             "for lib in dict.fromkeys((p['platlib'], p['purelib'])):\n"
             "    d = os.path.join(lib, 'lightsim2grid', 'share', 'cmake', 'lightsim2grid_core')\n"
             "    if os.path.isdir(d):\n"
             "        print(d)\n"
             "        break\n"],
            cwd=tmp, capture_output=True, text=True)
    found = fallback.stdout.strip()
    if fallback.returncode == 0 and found:
        print(f"WARNING: 'import lightsim2grid' failed with '{python}', using the CMake package found in "
              f"its site-packages. The import error was:\n{proc.stderr.strip().splitlines()[-1]}",
              file=sys.stderr)
        return Path(found.splitlines()[-1])
    raise RuntimeError(
        f"Could not get the lightsim2grid CMake directory with '{python}':\n{proc.stderr.strip()}\n"
        "Install lightsim2grid in that interpreter, pass --python / --cmake-dir, "
        "or build against the headers of this checkout with --source-tree.")


def check_cmake_dir(cmake_dir):
    """Fail if the package is not there; warn if the installed core has no KLU."""
    cmake_dir = Path(cmake_dir)
    config = cmake_dir / "lightsim2grid_coreConfig.cmake"
    if not config.is_file():
        raise RuntimeError(f"{config} not found: '{cmake_dir}' is not a lightsim2grid_core CMake package.")
    targets = cmake_dir / "lightsim2grid_coreTargets.cmake"
    has_klu = targets.is_file() and "KLU_SOLVER_AVAILABLE" in targets.read_text(errors="replace")
    if not has_klu:
        print("WARNING: KLU_SOLVER_AVAILABLE is not in lightsim2grid_coreTargets.cmake: the installed "
              "core was built without KLU, the *_KLU variants of the plugins will not register.",
              file=sys.stderr)
    return has_klu


def _run(cmd, cwd=None, env=None):
    print("+ " + " ".join(f'"{c}"' if " " in str(c) else str(c) for c in cmd), flush=True)
    return subprocess.run([str(c) for c in cmd], cwd=cwd, env=env).returncode


def _built_libraries(build_dir):
    return sorted(p for p in build_dir.rglob("*")
                  if p.is_file() and p.suffix in _LIB_SUFFIXES and "CMakeFiles" not in p.parts)


def build_plugins(plugins=None, *, cmake_dir=None, source_tree=False, python=None,
                  build_type="Release", jobs=None, clean=True, generator=None,
                  defines=None, test=False, keep_going=False):
    """Configure and build some example plugins (all of them by default).

    Parameters
    ----------
    plugins : list of str, optional
        Plugin folder names under examples/ (the "_algorithm" suffix may be left out).
    cmake_dir : str or Path, optional
        lightsim2grid_core CMake package to build against. Asked from the installed
        lightsim2grid when not given.
    source_tree : bool
        Build against the headers of this checkout instead of an installed package (the
        plugins' own fallback when no package is found).
    python : str, optional
        Interpreter whose lightsim2grid is used, both for ``cmake_dir`` and for ``test``.
        Defaults to the one running this script.
    build_type : str
        CMAKE_BUILD_TYPE, also the ``--config`` of multi-config generators (Visual Studio).
    jobs : int, optional
        Parallel build jobs (``cmake --build -j``). CMake's default when not given.
    clean : bool
        Remove ``<plugin>/build`` before configuring.
    generator : str, optional
        CMake generator (``-G``), e.g. "Ninja".
    defines : dict or list of "KEY=VALUE", optional
        Extra ``-D`` cache entries, passed to every plugin.
    test : bool
        Run each plugin's ``test_plugin.py`` after it is built.
    keep_going : bool
        Go on with the next plugin after a failure instead of stopping.

    Returns
    -------
    dict
        plugin name -> "ok" or the step that failed ("configure", "build", "test").
    """
    names = _resolve_names(plugins)
    python = python or sys.executable

    if source_tree:
        cmake_dir = None
        print("Building against the source tree of this checkout (no installed package).")
    else:
        cmake_dir = Path(cmake_dir) if cmake_dir else get_installed_cmake_dir(python)
        print(f"lightsim2grid CMake dir: {cmake_dir}")
        check_cmake_dir(cmake_dir)

    if isinstance(defines, dict):
        defines = [f"{k}={v}" for k, v in defines.items()]
    defines = list(defines or [])

    if shutil.which("cmake") is None:
        raise RuntimeError("cmake not found in PATH (pip install cmake, or install it system-wide).")

    status = {}
    for name in names:
        src = EXAMPLES_DIR / name
        build = src / "build"
        print(f"\n=== {name} ===", flush=True)
        if clean and build.exists():
            shutil.rmtree(build)

        cfg = ["cmake", "-S", src, "-B", build, f"-DCMAKE_BUILD_TYPE={build_type}"]
        if generator:
            cfg += ["-G", generator]
        if cmake_dir is not None:
            cfg.append(f"-DLIGHTSIM2GRID_CMAKE_DIR={cmake_dir}")
        cfg += [f"-D{d}" for d in defines]

        bld = ["cmake", "--build", build, "--config", build_type]
        if jobs:
            bld += ["-j", str(jobs)]

        step = None
        if _run(cfg) != 0:
            step = "configure"
        elif _run(bld) != 0:
            step = "build"
        else:
            for lib in _built_libraries(build):
                print(f"built: {lib}")
            test_file = src / "test_plugin.py"
            if test and test_file.is_file() and _run([python, test_file], cwd=EXAMPLES_DIR.parent) != 0:
                step = "test"

        status[name] = step or "ok"
        if step and not keep_going:
            break

    print("\n=== summary ===")
    for name in names:
        print(f"  {name:<30} {status.get(name, 'not run')}")
    return status


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Configure and build the lightsim2grid example plugins.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"Available plugins: {', '.join(available_plugins())}")
    parser.add_argument("plugins", nargs="*",
                        help="plugins to build (default: all); the '_algorithm' suffix may be left out")
    parser.add_argument("--list", action="store_true", help="list the plugins and exit")
    where = parser.add_mutually_exclusive_group()
    where.add_argument("--cmake-dir", help="lightsim2grid_core CMake package dir "
                                           "(default: lightsim2grid.get_cmake_dir() of the installed package)")
    where.add_argument("--source-tree", action="store_true",
                       help="build against the headers of this checkout, not an installed package")
    parser.add_argument("--python", help="interpreter whose lightsim2grid is used (default: this one)")
    parser.add_argument("--build-type", default="Release", help="CMAKE_BUILD_TYPE / --config (default: Release)")
    parser.add_argument("-j", "--jobs", type=int, help="parallel build jobs")
    parser.add_argument("--no-clean", dest="clean", action="store_false",
                        help="keep the existing build folders (incremental build)")
    parser.add_argument("-G", "--generator", help="CMake generator, e.g. Ninja")
    parser.add_argument("-D", dest="defines", action="append", default=[], metavar="KEY=VALUE",
                        help="extra CMake cache entry, passed to every plugin (repeatable)")
    parser.add_argument("--test", action="store_true", help="run each plugin's test_plugin.py after building it")
    parser.add_argument("-k", "--keep-going", action="store_true", help="do not stop at the first failure")
    args = parser.parse_args(argv)

    if args.list:
        print("\n".join(available_plugins()))
        return 0
    try:
        status = build_plugins(args.plugins, cmake_dir=args.cmake_dir, source_tree=args.source_tree,
                               python=args.python, build_type=args.build_type, jobs=args.jobs,
                               clean=args.clean, generator=args.generator, defines=args.defines,
                               test=args.test, keep_going=args.keep_going)
    except (RuntimeError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0 if all(s == "ok" for s in status.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
