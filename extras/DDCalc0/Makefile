############################################################
# Makefile for LUX rate calculation and likelihood         #
# routines.                                                #
#                                                          #
# Simple usage:                                            #
#   make all                                               #
#                                                          #
#   C Savage     University of Utah     2014               #
#                Nordita                2014               #
#                                                          #
############################################################

#################### INITIALIZATION ########################

# Generate variables needed for substitution.
empty :=
space := ${empty} ${empty}
lpar  := ${empty}(${empty}
rpar  := ${empty})${empty}


#################### VERSION ###############################

# Software release version.  Extracted from main code file.
DDCALC0_VERSION=$(shell grep "VERSION_STRING =" DDCalc0.f90 | cut -d \' -f2)
#$(info DDCalc0 version is $(DDCALC0_VERSION))


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
# for those compilers.  Any passed FC really only needs
# 'ifort' or 'gfortran' in the name, e.g. gfortran-4.9
# would be valid.
#ifeq (,$(filter $(FC),ifort gfortran))
ifeq (,$(findstring ifort,$(FC))$(findstring gfortran,$(FC)))
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
  ifneq (,$(findstring ifort,$(FC)))
    FFLAGS := -fast
  else ifneq (,$(findstring gfortran,$(FC)))
    FFLAGS := -O3 -fno-range-check
  else
    $(error Invalid compiler: $(FC).  Must use ifort or gfortran.)
  endif
endif

# Ensure required compiler flags
ifneq (,$(findstring gfortran,$(FC)))
  # disable compile-time range checking due to large integer
  # constants in math routines (allows for INTEGER*8)
  ifeq (,$(findstring -fno-range-check,$(FFLAGS)))
    FFLAGS += -fno-range-check
  endif
endif

# Flags used for shared object library.
# Shared library code must be compiled with -fPIC, which is
# a valid flag for both ifort and gfortran.  Not explicitly
# added to FFLAGS as there may be performance penalties in
# binaries (where it is not necessary), mainly arising from
# a reduction in compilation optimizations.
SOFFLAGS := $(FFLAGS)
ifeq (,$(findstring -fPIC,$(FFLAGS)))
  SOFFLAGS += -fPIC
endif


#################### FILES / TARGETS #######################

#f90objects := DDCalc0.o
#f90sources := $(f90objects:.o=.f90)
#f90mods := $(shell echo $(f90objects:.o=.mod) | tr A-Z a-z)

# Main programs
fprograms := DDCalc0run

# Module files (lower case as that is usual output filename)
fmodules := ddcalc0.mod

# Objects to be compiled
# e.g. file.o
fobjects := DDCalc0.o

# Libraries
statlib   := libDDCalc0.a
sharedlib := libDDCalc0.so
libraries := $(statlib) $(sharedlib)

# Include files
# e.g. file.h
fincludes := 

# Source files
fsources := $(fincludes) $(fobjects:.o=.f90) $(fprograms:=.f90)

# Additional files to include in tar file
extrafiles := Makefile DDCalc0run.use README.fortran DDCalc0.hpp

# Tar file
tarfile := DDCalc0.tar.gz

# Distribution tar file and subdirectory
distfile := DDCalc0-$(DDCALC0_VERSION).tar.gz
DIST_DIR := DDCalc0-$(DDCALC0_VERSION)

# Files to remove when cleaning.
# Note some compilers leave .dSYM files/directories for debugging.
# All generated files except for PDF and tar files.
cleanfiles := $(fprograms) $(fprograms:=.dSYM) \
              $(fobjects) $(fobjects:.o=.dSYM) \
              $(fmodules) \
              $(libraries)
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
$(fprograms): DDCalc0.f90

# Specific dependencies
#program1: object1.o


#################### RULES #################################

# Ensure that "all" is the default target.
#.DEFAULT_GOAL := all

# Default: build only the programs
.DEFAULT_GOAL := bin

# Phony targets
.PHONY: all bin lib doc tar dist FORCE \
        clean cleaner cleanest cleanish \
        clean-tar clean-dist \
        ifort gfortran

# Default is to build programs
all: $(fobjects) $(libraries) $(fprograms)

# Binaries (programs)
bin: $(fprograms)

# Libraries
lib: $(libraries)

# Define a do-nothing rule for dummy targets
# (avoids "nothing to be done" messages)
ifort gfortran:
	@:

# Rule for building programs
$(fprograms): % : %.f90
	@echo "===== Compiling $@ ====="
#	$(FC) $(FFLAGS) -o $@ $< $(filter %.o,$^)
	$(FC) $(FFLAGS) $(filter-out $<,$^) -o $@ $<

# Rule for building objects
$(fobjects): %.o : %.f90
	@echo "===== Compiling $@ ====="
	$(FC) $(FFLAGS) -o $@ -c $<

# Rule for building static library
$(statlib): $(fobjects)
	@echo "===== Generating $@ ====="
	ar rcs $@ $^

# Rule for building shared library
# We recompile from source to ensure -fPIC flag
# regardless if that flag is used elsewhere.
#$(sharedlib): $(fobjects)
#	@echo "===== Generating $@ ====="
#	$(FC) -shared -o $@ $^
$(sharedlib): $(fobjects:.o=.f90)
	@echo "===== Generating $@ ====="
	$(FC) $(SOFFLAGS) -shared -o $@ $^

# Generate tar archive
tar: $(tarfile)
$(tarfile): $(fsources) $(extrafiles)
	@echo "===== Generating $@ ====="
#	tar -czvf $@ $^
	tar -czf $@ $^

# Generate archive file
dist: $(distfile)
$(distfile): $(fsources) $(extrafiles)
	@echo "===== Generating $@ ====="
	mkdir -p ".archive-temp/$(DIST_DIR)"
	cp -al --parents $^ ".archive-temp/$(DIST_DIR)/"
#	tar -czvf $@ -C ".archive-temp" "$(DIST_DIR)"
	tar -czf $@ -C ".archive-temp" "$(DIST_DIR)"
	rm -rf ".archive-temp"

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

# Remove distributable tar file
clean-dist: FORCE
	@echo "===== Deleting distributable files ====="
	-rm -f $(distfile)

# Clean targets
# Files/directories to remove are explicitly given already,
# so we remove dependencies here.
clean:
cleaner: #clean
cleanest: #cleaner clean-tar
cleanish:


#################### END ###################################
