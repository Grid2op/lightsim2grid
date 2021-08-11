# Information
First, and most importantly, this project has been largely inspired from: 
https://github.com/jlblancoc/suitesparse-metis-for-windows

The goal of the file in this subdirectory is to be able to build suitesparse with cmake. 
Doing this allows to use the faster KLU solver even on windows based platform.

# How to build with cmake ?

## Requirements
This only works if you have a working `cmake` on your system. `cmake` is a free
software that you can download and install.

Of course, you also need to have a working compiler.

## Building SuiteSparse
You need to `cd` into the correct directory, which is the `build_cmake` subdirectory 
of this repo (`lightsim2grid/build_cmake`). Then:

1) you need to make sure to have clone the submodule, and in particular the
  "SuiteSparse" submodule, otherwise this will not work
2) generate the right sources (small hack to get cmake to work correctly):
   `python3 generate_c_files.py`
3) create a "build" subdirectory, and cd there: `mkdir build; cd build`
4) prepare the compilation with cmake: (on windows) 
  `cmake -DCMAKE_INSTALL_PREFIX=..\built -DCMAKE_BUILD_TYPE=Release ..` 
   **NB** It is mandatory to install the libraries there, otherwise they will not be detected
   by the installation script of lightsim2grid
5) compile the SuiteSparse package: `cmake --build . --config Release`
6) "install" it (which means: "move the libraries and the header at the right place): 
   `cmake --build . --config Release --target install`
   
**NB** Some of these files relies on symbolic links "symlink" and you might need root privilege on windows to
build them. We will try to have a workaround for that at some point.

**NB** On windows, the anti virus program (*eg* Avast, Windows Defender, MacAffee, etc.) might cause trouble during
the compilation. You might need to disable it temporarily during the compilation, or to add the proper exception 
to visual c++, cmake and possibly other programs.
