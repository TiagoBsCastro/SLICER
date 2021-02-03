# executable name
PREFIXDIR = /u/$(USER)
PROG = $(PREFIXDIR)/SLICER/bin/SLICER.EXE
OPTFLAGS = -DFixedPLCVertex # This optional flag forces the first snapshot to be centered around (0.5, 0.5, 0.5)

MAIN = slicer-v2.cpp utilities.cpp cosmology.cpp densitymaps.cpp gadget2io.cpp writeplc.cpp data.cpp

LDFLAGS += -Wl,-rpath -lstdc++ -lgsl -lgslcblas -lm -lCCfits -lcfitsio

ALLFLAGS =  -I$(PREFIXDIR)/include/ \
            -I$(PREFIXDIR)/include/gsl/ \
            -I$(PREFIXDIR)/include/CCfits/ \
            -I/opt/cluster/openmpi/3.1.6/gnu/9.3.0/include/ \
            -I./

LIBS = -L$(PREFIXDIR)/lib/ \
       -L/opt/cluster/openmpi/3.1.6/gnu/9.3.0/lib/

DEBUG = -g -O3 -fast -Wpedantic
CC = mpic++
RM = rm -f -r
OBJ = $(SOURCES:.cpp=.o)

CFLAGS=-O3 -g -fPIC -std=c++11

default: main
main:
	$(CC) $(CFLAGS) ${ALLFLAGS} $(MAIN) ${LIBS} ${LDFLAGS} ${OPTFLAGS} -o ${PROG}
clean:
	$(RM) $(PROG) $(OBJ) *~
