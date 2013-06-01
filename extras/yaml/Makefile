MAIN=example

CC = g++
LN = $(CC)
LNFLAGS = -lyaml-cpp
CFLAGS = -O3 -Wall -g

all: $(MAIN)

%.o: %.f
	$(CC) $(CFLAGS) -c $< -o $@

$(MAIN): $(MAIN).o
	$(LN) $(CFLAGS) $< -o $@ $(LNFLAGS)
