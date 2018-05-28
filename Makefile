# executable name
PROG = bin/MapSim-v5.10.4.Pinocchio

MAIN = main-v5.10.4.cpp  util.cpp

LDFLAGS += -Wl,-rpath  -lstdc++ -lgsl -lgslcblas -lm -lCCfits -lcfitsio

ALLFLAGS =  -I/$(SCRATCH)/include/ \
            -I/$(SCRATCH)/include/gsl/ \
            -I/$(SCRATCH)/include/CCfits/ \
	    -I./

LIBS = -L$(SCRATCH)/lib/ 

DEBUG = -g -O3 -fast -Wpedantic
CC = g++
RM = rm -f -r
OBJ = $(SOURCES:.cpp=.o)

CFLAGS=-O3 -g -fPIC -std=c++11

default: main
main:
	$(CC) $(CFLAGS) ${ALLFLAGS} $(MAIN) ${LIBS} ${LDFLAGS} -o ${PROG}
clean:
	$(RM) $(PROG) $(OBJ) *~
