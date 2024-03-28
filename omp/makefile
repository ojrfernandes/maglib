FIO_DIR=/home/jfernandes/Software/fusion-io/build

CXX = mpicxx
CFLAGS = -Wall -fPIC -O2 -fopenmp -g

LD = $(CXX)
LDFLAGS = -fPIC -Wl,-rpath,/opt/local/lib -fopenmp -g

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

fio_omp: fio_omp.cpp
	$(CXX) fio_omp.cpp -o fio_omp $(CFLAGS) $(INCLUDE) $(LDFLAGS) $(LIBS)

fio_omp_02 : fio_omp_02.o
	$(LD) $(LDFLAGS) $< $(LIBS) -o $@

fio_omp_02.o : fio_omp_02.cpp
	$(CXX) $< -c $(CFLAGS) $(INCLUDE) -o $@

korc_example : korc_example.o
	$(LD) $(LDFLAGS) $< $(LIBS) -o $@

korc_example.o : korc_example.cpp
	$(CXX) $< -c $(CFLAGS) $(INCLUDE) -o $@

clean :
	rm -f *.o *~ fio_omp fio_omp_02 korc_example