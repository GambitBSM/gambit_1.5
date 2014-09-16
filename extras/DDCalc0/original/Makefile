############################################################
# Makefile for LUX rate calculation and likelihood         #
# routines.                                                #
#                                                          #
# Simple usage:                                            #
#   make all                                               #
#                                                          #
#   C Savage     University of Utah     2014               #
#                                                          #
############################################################

#################### INITIALIZATION ########################

# Generate variables needed for substitution.
empty :=
space := ${empty} ${empty}
lpar  := ${empty}(${empty}
rpar  := ${empty})${empty}


#################### COMPILER ##############################

# The fortran compiler can be explicitly set here, along
# with the compilation options.  If not explicitly set, make
# will search for ifort and gfortran, in that order, and use
# one of those.  The compiler can also be be chosen by
# adding a 'ifort' or 'gfortran' target to the make command.

# This section simply sets FC and FFLAGS, so if you have
# difficulty running this makefile, simply set those two
# variables and remove everything else.

# Intel fortran compiler
#FC := ifort
#FFLAGS := -fast

# GNU fortran compiler
#FC := gfortran
#FFLAGS := -O3 -fno-range-check

# Fortran compiler (search)
# If FC not already set to ifort or gfortran, search
# for those compilers.
ifeq (,$(filter $(FC),ifort gfortran))
  FC := $(or $(shell which ifort 2> /dev/null),\
             $(shell which gfortran 2> /dev/null),\
             $(error Could not find ifort or gfortran compiler \
                     in the path.  A valid compiler must be given.))
  FC := $(notdir $(FC))
endif

# Explicit compiler, specified by additional 'ifort' or
# 'gfortran' target.
ifneq (,$(findstring ifort,$(MAKECMDGOALS)))
  FC := ifort
endif
ifneq (,$(findstring gfortran,$(MAKECMDGOALS)))
  FC := gfortran
endif

# Fortran compiler flags
ifeq (,$(FFLAGS))
  ifeq ($(FC),ifort)
    FFLAGS := -fast
  else ifeq ($(FC),gfortran)
    # disable range checking due to large integer constants
    # in math routines (allows for INTEGER*8)
    FFLAGS := -O3 -fno-range-check
  else
    $(error Invalid compiler: $(FC).  Must use ifort or gfortran.)
  endif
endif


#################### FILES / TARGETS #######################

#f90objects := DDCalc0.o
#f90sources := $(f90objects:.o=.f90)
#f90mods := $(shell echo $(f90objects:.o=.mod) | tr A-Z a-z)

# Main programs
fprograms := DDCalc0

# Module files (lower case as that is usual output filename)
fmodules := ddcalc0.mod

# Objects to be compiled
# e.g. file.o
fobjects := 

# Include files
# e.g. file.h
fincludes := 

# Source files
fsources := $(fincludes) $(fobjects:.o=.f90) $(fprograms:=.f90)

# Additional files to include in tar file
extrafiles := Makefile DDCalc0.use README.fortran

# Tar file
tarfile := DDCalc0.tar.gz

# Files to remove when cleaning.
# Note some compilers leave .dSYM files/directories for debugging.
# All generated files except for PDF and tar files.
cleanfiles := $(fprograms) $(fprograms:=.dSYM) \
              $(fobjects) $(fobjects:.o=.dSYM) \
              $(fmodules)
# ...also tar file.
cleanerfiles := $(cleanfiles) $(tarfile)
# ...and everything else.
cleanestfiles := $(cleanerfiles)
# Remove only build files that are not necessary for using program.
cleanishfiles := $(fprograms:=.dSYM) $(fobjects:.o=.dSYM)


#################### DEPENDENCIES ##########################

# General include file dependencies
#$(fobjects): file.h

# Some objects and programs have additional include files
#file.o: file_com.h

# Program dependencies
#$(fprograms): object1.o

# Specific dependencies
#program1: object1.o


#################### RULES #################################

# Ensure that "all" is the default target.
.DEFAULT_GOAL := all

# Phony targets
.PHONY: all doc tar FORCE \
        clean cleaner cleanest cleanish \
        clean-tar   \
        ifort gfortran

# Default is to build programs
all: $(fprograms)

# Define a do-nothing rule for dummy targets
# (avoids "nothing to be done" messages)
ifort gfortran:
	@:

# Rule for building programs
$(fprograms): % : %.f90
	@echo "===== Compiling $@ ====="
	$(FC) $(FFLAGS) -o $@ $< $(filter %.o,$^)

# Rule for building objects
$(fobjects): %.o : %.f90
	@echo "===== Compiling $@ ====="
	$(FC) $(FFLAGS) -o $@ -c $<

# Generate tar archive
tar: $(tarfile)
$(tarfile): $(fsources) $(extrafiles)
	@echo "===== Generating $@ ====="
#	tar -czvf $@ $^
	tar -czf $@ $^

# Remove files
clean: FORCE
	@echo "===== Cleaning ====="
	-rm -rf $(cleanfiles)
cleaner: FORCE
	@echo "===== Cleaning ====="
	-rm -rf $(cleanerfiles)
cleanest: FORCE
	@echo "===== Cleaning ====="
	-rm -rf $(cleanestfiles)
cleanish: FORCE
	@echo "===== Cleaning ====="
	-rm -rf $(cleanishfiles)

# Remove tar file
clean-tar: FORCE
	@echo "===== Deleting tar file ====="
	-rm -f $(tarfile)

# Clean targets
# Files/directories to remove are explicitly given already,
# so we remove dependencies here.
clean:
cleaner: #clean
cleanest: #cleaner clean-tar
cleanish:


#################### END ###################################
