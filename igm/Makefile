# executable name 
PROG = ./executable_igm

MAIN = main_igm.cpp readTNGParticle.cpp SED.cpp IGM_Ino14.cpp functions.cpp

LDFLAGS +=-Wl,-rpath -larmadillo -llapack -lblas -lstdc++


LIBS = -L/usr/lib/x86_64-linux-gnu/  -lgsl -lgslcblas  \
       -L/usr/lib/x86_64-linux-gnu/ \
       -L/usr/lib/x86_64-linux-gnu/ -lcfitsio -lCCfits \
       -L/usr/lib/x86_64-linux-gnu/ -lfftw3 -lm \
       -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lhdf5_cpp

# gsl, cfitsio, CCfits, fftw  
ALLFLAGS = -I/usr/local/include/gsl/ \
           -I/usr/include/ \
           -I/usr/local/include/CCfits/ \
           -I/usr/lib/x86_64-linux-gnu/ \
           -I/usr/local/include/ \
           -I//usr/include/eigen3/ \
           -I//usr/lib/x86_64-linux-gnu/hdf5/serial/ \
	   -I/usr/include/hdf5/serial/ \
           -I./

# 
DEBUG = -g -O2 

CC = /usr/bin/cc -std=c++11 -lstdc++ 

#
RM = rm -f -r
#
OBJ = $(SOURCES:.cpp=.o)
#

CFLAGS=-O2 -g -fPIC -DPIC

default: main
main: 
	$(CC) $(CFLAGS) ${ALLFLAGS} $(MAIN) ${LIBS} ${LDFLAGS} -o ${PROG}
clean:
	$(RM) $(PROG) $(OBJ) *~

