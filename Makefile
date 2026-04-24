# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of PyKLU2Grid, PyKLU2Grid a implements a c++ backend targeting the Grid2Op platform.


##################################
# this is the compilation of KLU
###################################
default: all

# CC = gcc
# CXX = g++
CC := $(if $(CC),$(CC),gcc)
CXX := $(if $(CXX),$(CXX),g++)

LIBPATH_REL = ./SuiteSparse
LIBPATH = $(realpath $(LIBPATH_REL))
PATHOUT = old

INCLUDE = -I$(LIBPATH)/SuiteSparse_config \
          -I$(LIBPATH)/CXSparse/Include \
          -I$(LIBPATH)/AMD/Include \
          -I$(LIBPATH)/BTF/Include \
          -I$(LIBPATH)/COLAMD/Include \
          -I$(LIBPATH)/KLU/Include

LIB = $(LIBPATH)/KLU/Lib/libklu.a \
      $(LIBPATH)/BTF/Lib/libbtf.a \
      $(LIBPATH)/AMD/Lib/libamd.a \
      $(LIBPATH)/COLAMD/Lib/libcolamd.a \
      $(LIBPATH)/CXSparse/Lib/libcxsparse.a \
      $(LIBPATH)/SuiteSparse_config/libsuitesparseconfig.a

clean:
	(cd $(LIBPATH)/SuiteSparse_config/ && make clean)
	(cd $(LIBPATH)/CXSparse/ && make clean)
	(cd $(LIBPATH)/AMD/ && make clean)
	(cd $(LIBPATH)/BTF/ && make clean)
	(cd $(LIBPATH)/COLAMD/ && make clean)
	(cd $(LIBPATH)/KLU/ && make clean)

purge:
	(cd $(LIBPATH)/SuiteSparse_config/ && make purge)
	(cd $(LIBPATH)/CXSparse/ && make purge)
	(cd $(LIBPATH)/AMD/ && make purge)
	(cd $(LIBPATH)/BTF/ && make purge)
	(cd $(LIBPATH)/COLAMD/ && make purge)
	(cd $(LIBPATH)/KLU/ && make purge)

distclean: purge

all: prelude

SS_CMAKE_OPTIONS = -DSUITESPARSE_REQUIRE_BLAS=OFF -DSUITESPARSE_USE_OPENMP=OFF -DBLA_VENDOR=Generic -DCMAKE_POSITION_INDEPENDENT_CODE=ON

# SuiteSparse v7.x packages no longer have a Lib/ sub-Makefile; their top-level
# Makefile delegates to cmake (building into <pkg>/build/).  After each build we
# copy the resulting static library into the legacy <pkg>/Lib/ location so that
# the lightsim2grid CMake "strategy 3" detector (which looks for <pkg>/Lib/lib*.a)
# continues to work.
prelude:
	(cd $(LIBPATH)/SuiteSparse_config/ && CC=$(CC) make CMAKE_OPTIONS="$(SS_CMAKE_OPTIONS)")
	cp $(LIBPATH)/SuiteSparse_config/build/libsuitesparseconfig.a $(LIBPATH)/SuiteSparse_config/libsuitesparseconfig.a
	(cd $(LIBPATH)/CXSparse/ && CC=$(CC) make CMAKE_OPTIONS="$(SS_CMAKE_OPTIONS)")
	mkdir -p $(LIBPATH)/CXSparse/Lib
	find $(LIBPATH)/CXSparse/build -maxdepth 2 -name "libcxsparse*.a" 2>/dev/null | head -1 | xargs -I{} cp {} $(LIBPATH)/CXSparse/Lib/libcxsparse.a
	(cd $(LIBPATH)/AMD/ && CC=$(CC) make CMAKE_OPTIONS="$(SS_CMAKE_OPTIONS)")
	mkdir -p $(LIBPATH)/AMD/Lib
	find $(LIBPATH)/AMD/build -maxdepth 2 -name "libamd*.a" 2>/dev/null | head -1 | xargs -I{} cp {} $(LIBPATH)/AMD/Lib/libamd.a
	(cd $(LIBPATH)/BTF/ && CC=$(CC) make CMAKE_OPTIONS="$(SS_CMAKE_OPTIONS)")
	mkdir -p $(LIBPATH)/BTF/Lib
	find $(LIBPATH)/BTF/build -maxdepth 2 -name "libbtf*.a" 2>/dev/null | head -1 | xargs -I{} cp {} $(LIBPATH)/BTF/Lib/libbtf.a
	(cd $(LIBPATH)/COLAMD/ && CC=$(CC) make CMAKE_OPTIONS="$(SS_CMAKE_OPTIONS)")
	mkdir -p $(LIBPATH)/COLAMD/Lib
	find $(LIBPATH)/COLAMD/build -maxdepth 2 -name "libcolamd*.a" 2>/dev/null | head -1 | xargs -I{} cp {} $(LIBPATH)/COLAMD/Lib/libcolamd.a
	(cd $(LIBPATH)/KLU/ && CC=$(CC) make CMAKE_OPTIONS="$(SS_CMAKE_OPTIONS)")
	mkdir -p $(LIBPATH)/KLU/Lib
	find $(LIBPATH)/KLU/build -maxdepth 2 -name "libklu*.a" ! -name "*cholmod*" 2>/dev/null | head -1 | xargs -I{} cp {} $(LIBPATH)/KLU/Lib/libklu.a

##################################
# this is the documentation
#################################
# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = docs
BUILDDIR      = documentation

doc_help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)


html: Makefile
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# cmd windows: > sphinx-build -M html docs documentation
# cmd for pdf > make latexpdf