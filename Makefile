# executable name
PREFIXDIR = /home/$(USER)
PROG = $(PREFIXDIR)/bin/SLICER.EXE

MAIN = slicer-v1.cpp util.cpp utilities.cpp cosmology.cpp densitymaps.cpp

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
