############################################################
# Minimal makefile for building the shared library         # 
############################################################

FC := gfortran

FFLAGS := -O3 -fno-range-check -shared -fPIC

library := DDCalc0.so

fsources := DDCalc0.f90

all: $(fsources)
	$(FC) $(FFLAGS) -o $(library) $(fsources)

clean:
	-rm -rf $(library)
