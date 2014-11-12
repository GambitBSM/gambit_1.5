CC = g++ 
FC = gfortran
CFLAGS = -O3 -Wall -Werror=uninitialized -g -fPIC -pedantic -fPIC -std=c++0x

BUILD = build
SRC = src
INCLUDE = include
LIB = lib
EXAMPLES = examples

CPPSOURCES:=$(wildcard $(SRC)/*.cpp)
CPPOBJECTS:=$(addprefix $(BUILD)/,$(notdir $(patsubst %.cpp,%.o,$(CPPSOURCES))))
CCSOURCES:=$(wildcard $(SRC)/*.cc)
CCOBJECTS:=$(addprefix $(BUILD)/,$(notdir $(patsubst %.cc,%.o,$(CCSOURCES))))
OBJECTS = $(CPPOBJECTS) $(CCOBJECTS)
TARGETS = pythia8_example #simple_example

# The following should really be done in a configure/cmake script rather than this makefile.
PY8_CXXFLAGS := $(shell pythia8-config --cxxflags --hepmc)
PY8_LDFLAGS := $(shell pythia8-config --libs --hepmc)
ROOT_CXXFLAGS := $(shell root-config --cflags)
ROOT_LDFLAGS := $(shell root-config --libs)
FJ_CXXFLAGS := $(shell fastjet-config --cxxflags)
FJ_LDFLAGS := $(shell fastjet-config --libs)
FJ_DIR := $(shell fastjet-config --prefix)
ifndef FJ_CXXFLAGS
  $(info ***!!!***)
  $(info Fastjet compile flags failed to configure.)
  $(info Please include /where/fastjet/is/bin in your $$PATH variable.)
  $(info You may also need to include /where/fastjet/is/lib in your)
  $(info $$LD_LIBRARY_PATH and $$DYLD_LIBRARY_PATH variables)
  $(info Cannot compile FastSim without these flags set!)
else
  $(info )
  $(info -- Fastjet Flags --)
  $(info Compile > $(FJ_CXXFLAGS))
  $(info Link    > $(FJ_LDFLAGS))
  $(info )
endif

.PHONY: clean $(TARGETS)

default: pythia8_example

all: libfastsim.so $(TARGETS)

$(BUILD)/%.o : $(SRC)/%.cpp
	$(CC) -c $(CFLAGS) $(FJ_CXXFLAGS) -I$(INCLUDE) $< -o $@

$(BUILD)/%.o : $(SRC)/%.cc
	$(CC) -c $(CFLAGS) $(FJ_CXXFLAGS) -I$(INCLUDE) $< -o $@

$(BUILD)/%.o : $(EXAMPLES)/%.cpp
	$(CC) -c $(CFLAGS) $(FJ_CXXFLAGS) $(PY8_CXXFLAGS) $(ROOT_CXXFLAGS) -I$(INCLUDE) $< -o $@

$(BUILD)/%.o : $(EXAMPLES)/%.cc
	$(CC) -c $(CFLAGS) $(FJ_CXXFLAGS) $(PY8_CXXFLAGS) $(ROOT_CXXFLAGS) -I$(INCLUDE) $< -o $@

libfastsim.so: $(OBJECTS)
	$(CC) -shared -o $(LIB)/$@ $^ $(FJ_LDFLAGS) -lgsl -lgslcblas

pythia8_pgun:  libfastsim.so $(BUILD)/example-fastsim-py8-pgun.o 
	$(CC) $(CFLAGS) -I$(INCLUDE) $(PY8_LDFLAGS) $(ROOT_LDFLAGS) -L$(LIB) -lfastsim $(BUILD)/example-fastsim-py8-pgun.o -o $@

pythia8_example: libfastsim.so $(BUILD)/example-fastsim-py8.o 
	$(CC) $(CFLAGS) -I$(INCLUDE) $(PY8_LDFLAGS) $(ROOT_LDFLAGS) -L$(LIB) -lfastsim $(BUILD)/example-fastsim-py8.o -o $@

#simple_example: libfastsim.so $(BUILD)/fastsim_simple.o 
#	$(CC) $(CFLAGS) -I$(INCLUDE) $(PY8_LDFLAGS) $(ROOT_LDFLAGS) -L$(LIB) -lfastsim $(BUILD)/fastsim_simple.o -o $@

clean:
	rm -rf $(TARGETS) $(LIB)/libfastsim.so $(BUILD)/*.o

