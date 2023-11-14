# executable name 
PROG = ./create_image

MAIN = create_image.cpp functions.cpp

LDFLAGS +=-Wl,-rpath -larmadillo -llapack -lblas -lstdc++ #-DSWAP #-DCOSMOLIB #,/scisoft/armadillo-10.1.0/



LIBS = -L/usr/lib/x86_64-linux-gnu/  -lgsl -lgslcblas  \
       -L/usr/lib/x86_64-linux-gnu/ \
       -L/usr/lib/x86_64-linux-gnu/ -lcfitsio -lCCfits \
       -L/usr/lib/x86_64-linux-gnu/ -lfftw3 -lm \
       -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lhdf5_cpp
##       -L/usr/lib64 -lstdc++ \
       -L/usr/lib/python3.6/ -lpython3.6m \

##-L/home/scisoft/armadillo-10.1.0/
##       -L/usr/local/hdf5/lib/libhdf5.so.103 \

       ##-L/usr/local/hdf5/lib/ -lhdf5_cpp -lhdf5   
       #-L$(COSMOLIB_DIR)/ -lCosmoLib \
       #-L$(NR_DIR)/ -lNR \

# gsl, cfitsio, CCfits, fftw  
ALLFLAGS = -I/usr/local/include/gsl/ \
           -I/usr/include/ \
           -I/usr/include/CCfits/ \
           -I/usr/lib/x86_64-linux-gnu/ \
           -I/usr/local/include/ \
           -I//usr/include/eigen3/ \
	   -I//usr/lib/x86_64-linux-gnu/hdf5/serial/ \
	   -I/usr/include/hdf5/serial/ \
           -I./
#           -I/home/scisoft/armadillo-10.1.0/include/armadillo_bits/ \
#           -I/home/scisoft/hdf5-1.10.7/hdf5/include/ \
#-I/usr/include/python3.6m/ \

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

