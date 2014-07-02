CC = g++ 
FC = gfortran
CCEXTRA= -pipe -pg -Wextra -g 
LNFLAGS += -g

PY8_CXXFLAGS := $(shell pythia8-config --cxxflags)
PY8_LDFLAGS := $(shell pythia8-config --libs)

ROOT_CXXFLAGS := $(shell root-config --cflags)
ROOT_LDFLAGS := $(shell root-config --libs)

FJ_CXXFLAGS := $(shell fastjet-config --cxxflags)
CFLAGS += -O3 -Wall -Werror=uninitialized -g -fPIC -pedantic -fopenmp  -std=c++11

BACKENDDIR = ../../../modules/Backends
#FASTSIMDIR = /import/sydpp2/atlas/Aldos_Data/gambit/hepforge2/extras/fast_sim
FASTSIMDIR = ../../fast_sim

$(FASTSIMDIR)/src/%.o : $(FASTSIMDIR)/src/%.cpp
	$(CC) -c $(CFLAGS) $(FJ_CXXFLAGS) -I$(FASTSIMDIR)/include -fPIC $< -o $@
$(FASTSIMDIR)/src/%.o : $(FASTSIMDIR)/src/%.cc
	$(CC) -c $(CFLAGS) $(FJ_CXXFLAGS) -I$(FASTSIMDIR)/include -fPIC $< -o $@

FJ_DIR := $(shell fastjet-config --prefix)

#$(BACKENDDIR)/lib/libfastsim.so: $(FASTSIMDIR)/src/FastSim.o $(FASTSIMDIR)/src/DetectorResponse.o $(FASTSIMDIR)/src/fastsim_interface.o $(FASTSIMDIR)/src/jsoncpp.o
#libfastsim.so: $(FASTSIMDIR)/src/FastSim_Reader.o $(FASTSIMDIR)/src/FastSim_Log.o $(FASTSIMDIR)/src/FastSim.o $(FASTSIMDIR)/src/DetectorResponse.o $(FASTSIMDIR)/src/fastsim_interface.o $(FASTSIMDIR)/src/jsoncpp.o


libfastsim.so: $(FASTSIMDIR)/src/FastSim_Log.o $(FASTSIMDIR)/src/DetectorResponse.o $(FASTSIMDIR)/src/jsoncpp.o $(FASTSIMDIR)/src/FastSim.o $(FASTSIMDIR)/src/FastSim_Reader.o 


# Intention is to pull symbols from FastJet into the fastsim backend (same should happen for direct pythia event filling... and any analyses that use jet structure)
# Currently breaks formost symbols because symbols can't be copied unless libfastjet.a was unusually built with the -fPIC flag
	$(CC) -shared -L$(FJ_DIR)/lib -o $@ $^ -lgsl -lgslcblas -Wl,-Bstatic -lfastjettools -Wl,-Bstatic -lfastjet -Wl,-Bdynamic
#	$(CC) -shared -o $@ $^ $(FJ_LDFLAGS)

example: example-fastsim-py8.cpp $(FASTSIMDIR)/src/libfastsim.so
	$(CC) $(CFLAGS) $(FJ_CXXFLAGS) -I$(FASTSIMDIR)/include $(PY8_CXXFLAGS) $(ROOT_CXXFLAGS) -fPIC $(PY8_LDFLAGS) $(ROOT_LDFLAGS) -L./ -lfastsim example-fastsim-py8.cpp -o $@

#simple: 
#	$(CC) $(CFLAGS) $(FJ_CXXFLAGS) -I$(FASTSIMDIR)/include $(PY8_CXXFLAGS) $(ROOT_CXXFLAGS) -fPIC $(PY8_LDFLAGS) $(ROOT_LDFLAGS) -L./ -lfastsim fastsim_simple.cpp -o $@



.PHONY: clean
clean:
			rm *.o example simple
