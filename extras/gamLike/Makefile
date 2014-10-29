BUILD=./
SRC=src
LIB=./

CXX=g++
AR=ar

CXXFLAGS=
CPPFLAGS=
ARFLAGS=rv

LDLIBS=-lgsl -lblas

SOURCES=gamLike.cc fermi_ini.cc fermi_like.cc

OBJ = $(patsubst %.cc,$(BUILD)/%.o,$(SOURCES)) 

all: $(LIB)/gamLike.a $(LIB)/gamLike.so

$(BUILD)/%.o: $(SRC)/%.cc src/gamLike.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(LIB)/gamLike.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $^

$(LIB)/gamLike.so: $(OBJ)
	$(CXX) -shared -rdynamic -o $@ $^ $(LDLIBS)

.PHONY: clean tarball
clean:
	rm -f $(BUILD)/*.o
	rm -f $(LIB)/*.o
	rm -f $(BUILD)/*.so
	rm -f $(LIB)/*.so
	rm -f $(BUILD)/*.a
	rm -f $(LIB)/*.a

tarball:
	tar -czvf gamLike.tar.gz $(SRC)/* data/like/*
