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

# The fortran compiler is set here, along with the
# compilation options.  If the compiler is not specified via
# an environmental variable or as an argument to the make
# invocation (FC=<...>), make will search for ifort and
# gfortran, in that order, and use one of those.

# This section simply sets FC and FFLAGS, so if you have
# difficulty running this makefile, simply set those two
# variables and remove everything else.

# Intel fortran compiler
#FC := ifort
#FFLAGS := -fast

# GNU fortran compiler
#FC := gfortran
#FFLAGS := -O3 -fno-range-check

# Command used to search for ifort or gfortran compiler.
# Gives an error if neither is found.
GETFC = $(or $(shell which ifort 2> /dev/null),\
             $(shell which gfortran 2> /dev/null),\
             $(error Could not find ifort or gfortran compiler \
                     in the path.  A valid compiler must be given.))

# If FC is not set, search for ifort or gfortran compilers,
# in that order.  If neither is found, make stops with an
# error message.
ifneq (,$(findstring "$(origin FC)","undefined" "default"))
  #FC := $(or $(shell which ifort 2> /dev/null),\
  #           $(shell which gfortran 2> /dev/null),\
  #           $(error Could not find ifort or gfortran compiler \
  #                   in the path.  A valid compiler must be given.))
  FC := $(GETFC)
  FC := $(notdir $(FC))
endif

# Check if FC is empty, which is only possible if it was
# explicitly set that way via argument or environmental
# variable.  Search for ifort and gfortran in that case.
ifeq (,$(FC))
  $(info WARNING: FC was specified but empty.  Searching for new compiler...)
  # If FC passed as an argument, must override
  ifeq ("$(origin FC)","command line")
    override FC := $(GETFC)
    override FC := $(notdir $(FC))
  else
    FC := $(GETFC)
    FC := $(notdir $(FC))
  endif
  $(info $(empty)  Found: $(FC))
endif

# Check if using ifort or gfortran and give a warning
# otherwise.  These are the only supported compilers.
ifeq (,$(findstring ifort,$(FC))$(findstring gfortran,$(FC)))
  $(info $(empty)------------------------------------------------------------)
  $(info $(empty)  WARNING:)
  $(info $(empty)  Only the ifort $(lpar)Intel$(rpar) and gfortran \
                 $(lpar)GCC$(rpar) compilers have)
  $(info $(empty)  been tested with DDCalc0.  Proceeding with unsupported)
  $(info $(empty)  compiler $(FC).)
  $(info $(empty)------------------------------------------------------------)
  #$(info $(empty))
endif

# Fortran compiler flags.  Update only if not already specified.
ifneq (,$(findstring "$(origin FFLAGS)","undefined" "default"))
  ifneq (,$(findstring ifort,$(FC)))
    FFLAGS := -fast
  else ifneq (,$(findstring gfortran,$(FC)))
    FFLAGS := -O3 -fno-range-check
  else
    $(info WARNING: No default compilation flags for compiler $(FC).)
  endif
endif

# Ensure required compiler flags.
ifneq (,$(findstring gfortran,$(FC)))
  # disable compile-time range checking due to large integer
  # constants in math routines (allows for INTEGER*8)
  ifeq (,$(findstring -fno-range-check,$(FFLAGS)))
    # If FFLAGS passed as an argument, must override
    ifeq ("$(origin FFLAGS)","command line")
      override FFLAGS += -fno-range-check
    else
      FFLAGS += -fno-range-check
    endif
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

# C++ compiler and flags.  Only used for test program.
# CXX_FLIBS are libraries necessary for linking fortran
# routines (compiler specific).
ifneq (,$(findstring "$(origin CXX)","undefined" "default"))
  ifneq (,$(findstring ifort,$(FC)))
    CXX := icc
  else ifneq (,$(findstring gfortran,$(FC)))
    CXX := g++
  endif
endif
ifneq (,$(findstring "$(origin CXXFLAGS)","undefined" "default"))
  ifneq (,$(findstring icc,$(CXX)))
    CXXFLAGS := 
  else ifneq (,$(findstring g++,$(CXX)))
    CXXFLAGS := 
  endif
endif
ifneq (,$(findstring icc,$(CXX)))
  CXX_FLIBS := -lifcore
else ifneq (,$(findstring g++,$(CXX)))
  CXX_FLIBS := -lgfortran
else
  $(info WARNING: Unsupported C++ compiler $(CXX).  Build may fail due to)
  $(info unknown library dependencies.)
  CXX_FLIBS := 
endif

# Debugging
#$(info FC      =$(FC) ($(origin FC)))
#$(info FFLAGS  =$(FFLAGS) ($(origin FFLAGS)))
#$(info SOFFLAGS=$(SOFFLAGS) ($(origin SOFFLAGS)))


#################### FILES / TARGETS #######################

#f90objects := DDCalc0.o
#f90sources := $(f90objects:.o=.f90)
#f90mods := $(shell echo $(f90objects:.o=.mod) | tr A-Z a-z)

# Main programs
fprograms := DDCalc0run

# Example/test programs
ftestprograms := DDCalc0_exampleF
ctestprograms := DDCalc0_exampleC

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
cincludes := DDCalc0.hpp

# Source files
fsources := $(fincludes) $(fobjects:.o=.f90) $(fprograms:=.f90) $(ftestprograms:=.f90)
csources := $(cincludes) $(ctestprograms:=.cpp)
sources  := $(fsources) $(csources)

# Additional files to include in tar file
extrafiles := Makefile DDCalc0run.use \
              README LICENSE CITE CITE.bib \
              example_efficiencies.dat

# Tar file
tarfile := DDCalc0.tar.gz

# Distribution tar file and subdirectory
distfile := DDCalc0-$(DDCALC0_VERSION).tar.gz
DIST_DIR := DDCalc0-$(DDCALC0_VERSION)

# Files to remove when cleaning.
# Note some compilers leave .dSYM files/directories for debugging.
# All generated files except for PDF and tar files.
cleanfiles := $(fprograms) $(fprograms:=.dSYM) \
              $(ftestprograms) $(ftestprograms:=.dSYM) \
              $(ctestprograms) \
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
$(ftestprograms): DDCalc0.f90
$(ctestprograms): DDCalc0.o

# Specific dependencies
#program1: object1.o


#################### RULES #################################

# Ensure that "all" is the default target.
#.DEFAULT_GOAL := all

# Default: build only the programs
.DEFAULT_GOAL := bin

# Phony targets
.PHONY: all bin lib examples doc tar dist FORCE \
        clean cleaner cleanest cleanish \
        clean-tar clean-dist \
        ifort gfortran

# Default is to build programs
all: $(fobjects) lib bin examples
#all: $(fobjects) $(libraries) $(fprograms) \
#     $(ftestprograms) $(ctestprograms)

# Binaries (programs)
bin: $(fprograms)

# Libraries
lib: $(libraries)

# Example/test binaries (programs)
examples: $(ftestprograms) $(ctestprograms)

# Define a do-nothing rule for dummy targets
# (avoids "nothing to be done" messages)
ifort gfortran:
	@:

# Rule for building programs
$(fprograms) $(ftestprograms): % : %.f90
	@echo "===== Compiling $@ ====="
#	$(FC) $(FFLAGS) -o $@ $< $(filter %.o,$^)
	$(FC) $(FFLAGS) $(filter-out $<,$^) -o $@ $<

# Rule for building programs (C++)
$(ctestprograms): % : %.cpp
	@echo "===== Compiling $@ ====="
	$(CXX) $(CXXFLAGS) -o $@ $< $(filter %.o,$^) $(CXX_FLIBS)

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
$(tarfile): $(sources) $(extrafiles)
	@echo "===== Generating $@ ====="
#	tar -czvf $@ $^
	tar -czf $@ $^

# Generate archive file
dist: $(distfile)
$(distfile): $(sources) $(extrafiles)
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
