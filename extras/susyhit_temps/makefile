OBJS1 = bsg.o twoloophiggs.o suspect2.o
OBJS2 = sdecay.o lightst4bod.o Xvegas.o lightst3bod.o
OBJS3 = susylha.o hgaga.o hdecay.o

FF=gfortran

FC=$(FF)
FFLAGS= -O2 -fPIC

.f.o: 
	$(FC) $(FFLAGS) -c -finit-local-zero -fbacktrace -fno-align-commons $*.f

susyhit.so:$(OBJS1) $(OBJS2) $(OBJS3) 
	$(FC) -shared -o $@ $(OBJS1) $(OBJS2) $(OBJS3)

clean: 
	rm -f *.o
