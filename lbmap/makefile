SOURCE_DIR = /home/jfernandes/maglib/maglit/tcabr_lobes/src
OBJ_DIR = /home/jfernandes/Software/maglib/maglit/tcabr_lobes/obj

CXX = mpicxx
CFLAGS = -Wall -fPIC -O2

LD = $(CXX)
LDFLAGS = -fPIC -Wl,-rpath,/opt/local/lib

INCLUDE = -I/opt/local/include -I/usr/include
LIBS = -L/opt/local/lib/gcc9 -lgfortran -lstdc++

HDF5_LIBS = -L/home/jfernandes/CMake-hdf5-1.14.0/hdf5-1.14.0/hdf5/lib -lhdf5
HDF5_INCLUDE = -I/home/jfernandes/CMake-hdf5-1.14.0/hdf5-1.14.0/hdf5/include
NETCDF_LIBS = -lnetcdf

LAPACK = -llapack

LIBS := $(HDF5_LIBS) $(LIBS)

INCLUDE := $(HDF5_INCLUDE) $(INCLUDE)

# Source files
SRCS = $(SOURCE_DIR)/run.cpp $(SOURCE_DIR)/lobes.cpp $(SOURCE_DIR)/intersections.cpp

# List of object files
OBJS = $(OBJ_DIR)/run.o $(OBJ_DIR)/intersections.o $(OBJ_DIR)/lobes.o

# Default target
run: $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

# Compile source files into object files
$(OBJ_DIR)/%.o: $(SOURCE_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CFLAGS) $(INCLUDE) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJS) *~ run