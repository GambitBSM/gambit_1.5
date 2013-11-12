# Makefile for ge directory
# Author: R. Ruiz


.KEEP_STATE:

.SUFFIXES:      .o .f90


DINC =../../inc
MCMCINC=../../..

.f90.o:
	$(FC) $(FCFLAGS) -c -I$(DINC)  -I$(MCMCINC) $<

OBJ =  $(subst .f90,.o,  $(filter-out  $(special),$(wildcard *.f90) ))

all: $(OBJ)