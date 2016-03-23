# ============================================================================
# Name        : Makefile
# Author      : Ricardo Serrano
# Version     :
# Copyright   : This is CopyRighted to Ricardo Serrano
# Description : Makefile for Hello MPI World in Fortran
# ============================================================================

.PHONY: all clean

.SUFFIXES: .for

FC=gfortran -g3 -ggdb
LIBS=-framework vecLib -lpthread -lslatec
FLAGS=
FLAGS+=$(INCS) $(LIBS)

SOURCES=mtxutil.for gaussin.for matutil.for femutil.for uelutil.for asemutil.for postproc.for preproc.for apiro.for 

OBJECTS=$(SOURCES:.for=.o)

all: $(SOURCES) piro

piro: $(OBJECTS)
	$(FC) $(FLAGS) -o $@ $(OBJECTS)


.for.o:
	$(FC) $(INCS) -c $< -o $@

clean:
	rm -f piro *.mod *.o

echo:
	echo $(OBJS)
