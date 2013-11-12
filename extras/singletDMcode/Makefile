FC = ifort
FCFLAGS = -O3 -fpp -fPIC -parallel -openmp -lpthread -check bounds
export FC FCFLAGS
LIB_DM=dm/lib
LAPACKL = -llapack -lblas
OBJ = types.o numerical_routines.o jimlib.o dwarfs.o cta.o cmb.o driver.o
OBJ2 = types.o jimlib.o postprocess.o
OBJ3 = types.o jimlib.o postprocess_bfs.o
INCLUDE = .
PNAME1 = singlet
PNAME2 = postprocess
PNAME3 = postprocess_bfs
LIBS =  -L$(LIB_DM) -ldm $(LAPACKL)
TYPEDEF = typedef.o

default : all

all:  $(PNAME1) $(PNAME2) $(PNAME3)

%.o: %.f90 
	$(FC) $(FCFLAGS) -I$(INCLUDE) -c $*.f90

%.o: %.F90 
	$(FC) $(FCFLAGS) -I$(INCLUDE) -c $*.F90

$(PNAME1): $(TYPEDEF) dmlib $(OBJ) 
	$(FC) -o $(PNAME1) $(TYPEDEF) $(OBJ) $(LIBS) $(FCFLAGS)

$(PNAME2): $(OBJ2)  
	$(FC) -o $(PNAME2) $(OBJ2) $(FCFLAGS)

$(PNAME3): $(OBJ3)  
	$(FC) -o $(PNAME3) $(OBJ3) $(FCFLAGS)

dmlib: 
	cd dm ; $(MAKE) all


clean: 
	rm -f *.o *.mod *.d *.pc *.obj core $(PNAME1) $(PNAME2) $(PNAME3)

clean_dm:
	cd dm; $(MAKE) clean;

cleanall:
	$(MAKE) clean clean_dm

