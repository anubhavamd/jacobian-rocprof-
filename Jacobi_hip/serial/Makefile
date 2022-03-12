##**************************************************************************
##* Copyright (c) 2019, Advanced Micro Devices, Inc. All rights reserved.
##**************************************************************************

# Compilers
MPICC=$(PREP) mpic++
MPILD=$(PREP) mpic++

MPI_HOME=/usr/lib/openmpi/

# Flags
CFLAGS=-O3 -g -fPIC -std=c++11 -march=native -Wall
MPICFLAGS=-I${MPI_HOME}/include
MPILDFLAGS=$(shell mpic++ --showme:link)

# Description of binaries
BINDIR=.
JACOBI=$(BINDIR)/Jacobi
BINARIES=$(JACOBI)

DEPS=Jacobi.hpp \
		 defines.hpp

OBJS=JacobiSetup.o  \
		 JacobiRun.o    \
		 JacobiMain.o   \
		 HaloExchange.o \
		 Input.o        \
		 JacobiIteration.o \
		 Laplacian.o    \
		 Norm.o

# Commands
all: $(BINARIES)

HaloExchange.o: $(DEPS)  HaloExchange.cpp
	$(MPICC) $(MPICFLAGS) $(CFLAGS) -c HaloExchange.cpp -o HaloExchange.o

Input.o: $(DEPS)  Input.cpp
	$(MPICC) $(MPICFLAGS) $(CFLAGS) -c Input.cpp -o Input.o

JacobiMain.o: $(DEPS)  JacobiMain.cpp
	$(MPICC) $(MPICFLAGS) $(CFLAGS) -c JacobiMain.cpp -o JacobiMain.o

JacobiRun.o: $(DEPS)  JacobiRun.cpp
	$(MPICC) $(MPICFLAGS) $(CFLAGS) -c JacobiRun.cpp -o JacobiRun.o

JacobiSetup.o: $(DEPS)  JacobiSetup.cpp
	$(MPICC) $(MPICFLAGS) $(CFLAGS) -c JacobiSetup.cpp -o JacobiSetup.o

JacobiIteration.o: $(DEPS)  JacobiIteration.cpp
	$(MPICC) $(MPICFLAGS) $(CFLAGS) -c JacobiIteration.cpp -o JacobiIteration.o

Laplacian.o: $(DEPS)  Laplacian.cpp
	$(MPICC) $(MPICFLAGS) $(CFLAGS) -c Laplacian.cpp -o Laplacian.o

Norm.o: $(DEPS)  Norm.cpp
	$(MPICC) $(MPICFLAGS) $(CFLAGS) -c Norm.cpp -o Norm.o

$(JACOBI): $(OBJS)
	$(MPILD) $(MPICFLAGS) $(CFLAGS) -o $(JACOBI)  $(OBJS) $(MPILDFLAGS)

clean:
	rm -rf *.o *~ $(BINARIES)
