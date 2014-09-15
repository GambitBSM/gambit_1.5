############################################################
# Minimal makefile for building the shared library         # 
############################################################

FC = gfortran

FFLAGS = -O3 -fPIC

library = libDDCalc0.so

fobjects = DDCalc0.o

all: $(library)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(library): $(fobjects)
	$(FC) -shared -o $(library) $(fobjects)

clean:
	-rm -rf $(library) $(fobjects) *.mod
