# Makefile for ge directory
# Author: R. Ruiz


.KEEP_STATE:

.SUFFIXES:      .o .f

DINC =../../inc

.f.o:
	$(FC) $(FCFLAGS) -c -I$(DINC)  $<

OBJ =  $(subst .f,.o,  $(filter-out  $(special),$(wildcard *.f) ))

all: $(OBJ)