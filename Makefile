# executable name
PREFIXDIR = /data/tiago/
PROG = $(PREFIXDIR)/bin/MapSim-v5.10.6

MAIN = main-v5.10.6.cpp  util.cpp

LDFLAGS += -Wl,-rpath -lstdc++ -lgsl -lgslcblas -lm -lCCfits -lcfitsio -lcurl

ALLFLAGS =  -I/$(PREFIXDIR)/include/ \
            -I/$(PREFIXDIR)/include/gsl/ \
            -I/$(PREFIXDIR)/include/CCfits/ \
            -I/$(PREFIXDIR)/include/curl/ \
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
