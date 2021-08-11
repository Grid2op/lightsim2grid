# Information
First, and most importantly, this project has been largely inspired from: 
[suitesparse-metis-for-windows](https://github.com/jlblancoc/suitesparse-metis-for-windows) (we do not include the
compilation of "metis" because we don't need it for this project).

The goal of the files in this subdirectory is to be able to build some of the `SuiteSparse` modules with cmake 
(as opposed to the recommended compilation tool "make" used by this code). This allows to further speed up
`lightsim2grid` code with the use of the `KLU` linear solver (part of SuiteSparse) which is faster than the default
linear solver in `Eigen`.
 
Being able to compile SuiteSparse with cmake allows to benefit the faster `KLU` solver even on windows based platform, 
where "make" is not a viable option.

# How to build with cmake ?

## Requirements
This only works if you have a working `cmake` on your system. `cmake` is a free
software that you can download and install. We do not cover this (have a look at https://cmake.org/ 
for more information)

Of course, you also need to have a working compiler (for example microsoft visual studio on windows: 
https://visualstudio.microsoft.com/vs/features/cplusplus/ )

## Building SuiteSparse
You need to `cd` into the correct directory, which is the `build_cmake` subdirectory 
of this repo (`lightsim2grid/build_cmake`). Then:

1) you need to make sure to have clone the submodule, and in particular the
  "SuiteSparse" submodule, otherwise this will not work (`git submodule init; git submodule update`)
2) generate the right sources (small hack to get cmake to work correctly):
   `python3 generate_c_files.py` [might require admin rights for now...]
3) create a "build" subdirectory, and cd there: `mkdir build; cd build`
4) prepare the compilation with cmake: (on windows, adapt the command the other platforms) 
  `cmake -DCMAKE_INSTALL_PREFIX=..\built -DCMAKE_BUILD_TYPE=Release ..` 
   **NB** It is mandatory to install the libraries in `..\built`, otherwise they will not be detected
   by the installation script of `lightsim2grid` and this will be useless
5) compile the SuiteSparse package: `cmake --build . --config Release`
6) "install" it (which means: "move the libraries and the header at the right place): 
   `cmake --build . --config Release --target install`

## Install lightsim2grid

Now you can install lightsim2grid as usual (see the standard readme for more details):
1) `cd ..` (to go to the source tree of lightsim2grid)
2) `py -m pip install .`

## Test the installation is working

To assess whether this works or not you can now, from a python "shell" (really called "repl")

```python
from lightsim2grid.solver import KLUSolver
```

If you don't have any **errors** you are good to go.

**NB** you might encounter some **warnings** like 
"*Numba cannot be loaded. You will gain possibly massive speed if installing it*". To silence this warning, just install
the python package `numba`. This does not come from lightsim2grid, but is used by pandapower.
 
## Troubleshoot

Compiling on windows might be tricky. Here are some errors we encounter. Feel free to post a github issue if you
encounter some not listed there.

### Permission issue (read / write)
Some of these files relies on symbolic links "symlink" and you might need root privilege on windows to
build them. We will try to have a workaround for that at some point.

### Impossible to execute visual c++
On windows, the anti virus program (*eg* Avast, Windows Defender, MacAffee, etc.) might cause trouble during
the compilation. You might need to disable it temporarily during the compilation, or to add the proper exception 
to visual c++, cmake and possibly other programs.
