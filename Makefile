# Simple C Wrapper Makefile. Provide routine for solving linear system.

default: all

CC = gcc
CXX = g++

LIBPATH_REL = ./SuiteSparse
LIBPATH = $(realpath $(LIBPATH_REL))

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
	rm libpyklu.so

purge:
	(cd $(LIBPATH)/SuiteSparse_config/ && make purge)
	(cd $(LIBPATH)/CXSparse/Lib && make purge)
	(cd $(LIBPATH)/AMD/Lib/ && make purge)
	(cd $(LIBPATH)/BTF/Lib/ && make purge)
	(cd $(LIBPATH)/COLAMD/Lib/ && make purge)
	(cd $(LIBPATH)/KLU/Lib/ && make purge)
	rm libpyklu.so

distclean: purge

all: prelude pyKLU

prelude:
	(cd $(LIBPATH)/SuiteSparse_config/ && make CC=gcc)
	(cd $(LIBPATH)/CXSparse/Lib/ && make CC=gcc)
	(cd $(LIBPATH)/AMD/Lib/ && make CC=gcc)
	(cd $(LIBPATH)/BTF/Lib/ && make CC=gcc)
	(cd $(LIBPATH)/COLAMD/Lib/ && make CC=gcc)
	(cd $(LIBPATH)/KLU/Lib/ && make CC=gcc)

pyKLU: pyklu.c $(LIB)
	$(CC) $(INCLUDE) -fPIC -shared -o libpyklu.so pyklu.c $(LIB)
	# $(CXX) $(INCLUDE) -fPIC -shared -o libpyklu_cpp.so pyklu.cpp $(LIB)

