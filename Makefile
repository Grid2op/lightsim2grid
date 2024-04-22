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
	(cd $(LIBPATH)/CXSparse/Lib && make clean)
	(cd $(LIBPATH)/AMD/Lib/ && make clean)
	(cd $(LIBPATH)/BTF/Lib/ && make clean)
	(cd $(LIBPATH)/COLAMD/Lib/ && make clean)
	(cd $(LIBPATH)/KLU/Lib/ && make clean)

purge:
	(cd $(LIBPATH)/SuiteSparse_config/ && make purge)
	(cd $(LIBPATH)/CXSparse/Lib && make purge)
	(cd $(LIBPATH)/AMD/Lib/ && make purge)
	(cd $(LIBPATH)/BTF/Lib/ && make purge)
	(cd $(LIBPATH)/COLAMD/Lib/ && make purge)
	(cd $(LIBPATH)/KLU/Lib/ && make purge)

distclean: purge

all: prelude

prelude:
	(cd $(LIBPATH)/SuiteSparse_config/ && CC=$(CC) make)
	(cd $(LIBPATH)/CXSparse/Lib/ && CC=$(CC) make )
	(cd $(LIBPATH)/AMD/Lib/ && CC=$(CC) make)
	(cd $(LIBPATH)/BTF/Lib/ && CC=$(CC) make)
	(cd $(LIBPATH)/COLAMD/Lib/ && CC=$(CC) make)
	(cd $(LIBPATH)/KLU/Lib/ && CC=$(CC) make)

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