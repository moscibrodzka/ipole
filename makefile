#
# h5cc compiles for linking with HDF5 library
#
CC = h5cc -DH5_USE_16_API 
CFLAGS =  -fopenmp -I/usr/include -O3 -w
LDFLAGS = -lm -lgsl -lgslcblas 

SRCIPO = \
main.c image.c geodesics.c radiation.c tetrads.c ipolarray.c geometry.c \
model_tetrads.c model_radiation.c model_geodesics.c \
model_harm3d.c
#model_NT.c


OBJIPO = \
main.o image.o geodesics.o radiation.o tetrads.o ipolarray.o geometry.o \
model_tetrads.o model_radiation.o model_geodesics.o \
model_harm3d.o
#model_NT.o


ipole: $(OBJIPO) makefile 
	$(CC) $(CFLAGS) -o ipole $(OBJIPO) $(LDFLAGS)

$(OBJIPO) : makefile decs.h defs.h constants.h

clean:
	rm *.o 
cleanup:
	rm ipole*.ppm ipole.dat



