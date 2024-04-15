FIO_DIR=/home/jfernandes/Software/fusion-io/build

CXX = mpicxx
CFLAGS = -Wall -fPIC -O2

LD = $(CXX)
LDFLAGS = -fPIC -Wl,-rpath,/opt/local/lib

INCLUDE = -I/opt/local/include 
LIBS = -L/opt/local/lib/gcc9 -lgfortran -lstdc++ 

HDF5_LIBS = -L/home/jfernandes/Software/CMake-hdf5-1.14.0/hdf5-1.14.0/hdf5/lib -lhdf5
HDF5_INCLUDE = -I/home/jfernandes/Software/CMake-hdf5-1.14.0/hdf5-1.14.0/hdf5/include
NETCDF_LIBS = -lnetcdf

LAPACK = -llapack

LIBS := -L$(FIO_DIR)/lib -lfusionio -lm3dc1 \
		-Wl,-rpath,$(FIO_DIR)/lib \
		$(HDF5_LIBS) $(LIBS)

INCLUDE := -I$(FIO_DIR)/include $(HDF5_INCLUDE) $(INCLUDE)

test : test.o
	$(LD) $(LDFLAGS) $< $(LIBS) -o $@

test.o : test.cpp
	$(CXX) $< -c $(CFLAGS) $(INCLUDE) -o $@

clean :
	rm -f *.o *~ test