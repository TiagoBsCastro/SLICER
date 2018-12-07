# executable name
PREFIXDIR = /home/$(USER)
PROG = $(PREFIXDIR)/bin/MapSim-v7

#MAIN = main-v6.cpp  util.cpp cosmology.cpp utilities.cpp
MAIN = main-v7.cpp util_new.cpp utilities.cpp cosmology.cpp

LDFLAGS += -Wl,-rpath -lstdc++ -lgsl -lgslcblas -lm -lCCfits -lcfitsio

ALLFLAGS =  -I$(PREFIXDIR)/include/ \
            -I$(PREFIXDIR)/include/gsl/ \
            -I$(PREFIXDIR)/include/CCfits/ \
	          -I./

LIBS = -L$(PREFIXDIR)/lib/

DEBUG = -g -O3 -fast -Wpedantic
CC = mpic++
RM = rm -f -r
OBJ = $(SOURCES:.cpp=.o)

CFLAGS=-O3 -g -fPIC -std=c++11

default: main
main:
	$(CC) $(CFLAGS) ${ALLFLAGS} $(MAIN) ${LIBS} ${LDFLAGS} -o ${PROG}
clean:
	$(RM) $(PROG) $(OBJ) *~
