FIO_DIR = /home/jfernandes/Software/fusion-io/build
MAGLIT_DIR = /home/jfernandes/Software/maglib/maglit
SODE_DIR = /home/jfernandes/Software/maglib/sode
SOURCE_DIR = /home/jfernandes/Software/maglib/maglit/tcabr_footprints/source
OBJ_DIR = $(SOURCE_DIR)/obj

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
SRCS = $(SOURCE_DIR)/run.cpp $(MAGLIT_DIR)/maglit.cpp $(SODE_DIR)/sode.cpp $(SOURCE_DIR)/tcabr_collider.cpp $(SOURCE_DIR)/input_read.cpp $(SOURCE_DIR)/footprint.cpp

# List of object files
OBJS = $(OBJ_DIR)/run.o $(OBJ_DIR)/maglit.o $(OBJ_DIR)/sode.o $(OBJ_DIR)/tcabr_collider.o $(OBJ_DIR)/input_read.o $(OBJ_DIR)/footprint.o

# Default target
run: $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

# Compile source files into object files
$(OBJ_DIR)/%.o: $(SOURCE_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CFLAGS) $(INCLUDE) -c $< -o $@

# Compile maglit.cpp into object file
$(OBJ_DIR)/maglit.o: $(MAGLIT_DIR)/maglit.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CFLAGS) $(INCLUDE) -c $< -o $@

# Compile sode.cpp into object file
$(OBJ_DIR)/sode.o: $(SODE_DIR)/sode.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CFLAGS) $(INCLUDE) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJS) *~ run