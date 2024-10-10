FIO_DIR = /home/jfernandes/Software/fusion-io/build
MAGLIT_DIR = /home/jfernandes/Software/maglib/maglit
SODE_DIR = /home/jfernandes/Software/maglib/sode

CXX = mpicxx
CFLAGS = -Wall -fPIC -O2 -g -fopenmp 
#-fsanitize=thread

LD = $(CXX)
LDFLAGS = -fPIC -Wl,-rpath,/opt/local/lib -fopenmp

INCLUDE = -I/opt/local/include 
LIBS = -L/opt/local/lib/gcc9 -lgfortran -lstdc++ 

HDF5_LIBS = -L/home/jfernandes/Software/CMake-hdf5-1.14.0/hdf5-1.14.0/hdf5/lib -lhdf5
HDF5_INCLUDE = -I/home/jfernandes/Software/CMake-hdf5-1.14.0/hdf5-1.14.0/hdf5/include
NETCDF_LIBS = -lnetcdf

LAPACK = -llapack

LIBS := -L$(FIO_DIR)/lib -lfusionio -lm3dc1 \
		-Wl,-rpath,$(FIO_DIR)/lib \
		$(HDF5_LIBS) $(LIBS)

INCLUDE := -I$(FIO_DIR)/include $(HDF5_INCLUDE) $(INCLUDE) \
		-I$(MAGLIT_DIR) \
		-I$(SODE_DIR)

# Source files
SRCS = test.cpp $(MAGLIT_DIR)/maglit.cpp $(SODE_DIR)/sode.cpp ../tcabr_footprints/tcabr_collider.cpp

# List of object files
OBJS = $(SRCS:.cpp=.o)

# Default target
test: $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

# Compile source files into object files
$(MAGLIT_DIR)/%.o: $(MAGLIT_DIR)/%.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) -c $< -o $@

$(SODE_DIR)/%.o: $(SODE_DIR)/%.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJS) *~ test