# Optimization flags
ifeq ($(OPTIMIZE),1)
$(warning *** Generating optimized build! CUDA error checking is OFF! ***)
OPT_CFLAGS = -O2
OPT_NVCCFLAGS = 
else
# Debug mode is default
OPT_CFLAGS = -g
OPT_NVCCFLAGS = -DDEBUG 
endif

JSONCPP_SRC = contrib/jsoncpp-src-0.6.0-rc2/src/lib_json
JSONCPP_INC = contrib/jsoncpp-src-0.6.0-rc2/include

INCLUDE = -Isrc -I$(RATROOT)/include -I$(ROOTSYS)/include -I$(RATROOT)/src/stlplus -Icontrib/hemi -I/usr/local/cuda/include -I/opt/local/include -I/opt/cuda-5.0/include -I$(JSONCPP_INC)
CFLAGS = -DVERBOSE=true $(OPT_CFLAGS) $(INCLUDE)
GCCFLAGS = -Wall -Werror -Wno-unused-variable -ftrapv -fdiagnostics-show-option  # -Wunused-variable errors with HEMI macros
NOT_NVCC_CFLAGS =
NVCCFLAGS = -gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=\"sm_35,compute_35\" -use_fast_math $(OPT_NVCCFLAGS)
ROOTLIBS =  -lCore -lCint -lRIO -lMathCore -lHist -lGpad -lTree -lTree -lGraf -lm
LFLAGS = -L$(RATROOT)/lib -lRATEvent_$(RATSYSTEM) -L$(ROOTSYS)/lib $(ROOTLIBS) -L/opt/local/lib -lhdf5 -lhdf5_hl

# Mac hacks!
ARCH = $(shell uname)

ifndef CUDA_ROOT
$(warning *** CUDA_ROOT is not set, defaulting to CPU-only build ***)
GCC = g++ $(GCCFLAGS) $(NOT_NVCC_CFLAGS) -DHEMI_CUDA_DISABLE
CUDACC = $(CC)
CC = $(GCC)
else
GCC = g++ $(GCCFLAGS) $(NOT_NVCC_CFLAGS)
	ifeq ($(ARCH), Darwin)
		NVCCFLAGS := -m64 $(NVCCFLAGS)
		CUDA_LIBDIR = $(CUDA_ROOT)/lib
	else
		CUDA_LIBDIR = $(CUDA_ROOT)/lib64
	endif
CUDA_LFLAGS = -L$(CUDA_LIBDIR) -lcudart -lcurand
CC = $(CUDA_ROOT)/bin/nvcc $(NVCCFLAGS) -I$(CUDA_ROOT)/include
CUDACC = $(CC) -x cu
CC += --compiler-options "$(GCCFLAGS) -Wno-unused-function"
endif

OBJ_DIR = build
SOURCES = $(filter-out src/mcmc.cpp src/nll_kernels.cpp src/pdfz.cpp, $(wildcard src/*.cpp))
OBJECTS = $(SOURCES:src/%.cpp=$(OBJ_DIR)/%.o)
JSONCPP_SOURCES = $(wildcard $(JSONCPP_SRC)/*.cpp)
JSONCPP_OBJECTS = $(JSONCPP_SOURCES:$(JSONCPP_SRC)/%.cpp=$(OBJ_DIR)/jsoncpp/%.o)

# For unit test suit
SXMC_NO_MAIN_FUNCTION_OBJECTS = $(filter-out build/sensitivity.o, $(OBJECTS) $(JSONCPP_OBJECTS) build/mcmc.o build/nll_kernels.o build/pdfz.o)
TEST_SOURCES = $(wildcard test/*.cpp)
TEST_OBJECTS = $(TEST_SOURCES:test/%.cpp=$(OBJ_DIR)/test/%.o)

EXE = bin/sensitivity

ifndef RATROOT
$(error RATROOT is not set)
endif

ifndef ROOTSYS
$(error ROOTSYS is not set)
endif

all: $(OBJ_DIR)/mcmc.o $(OBJ_DIR)/nll_kernels.o $(OBJ_DIR)/pdfz.o $(OBJECTS) $(JSONCPP_OBJECTS) $(EXE)

.PHONY: doc test

clean:
	-$(RM) build/*.o build/test/*.o build/jsoncpp/*.o bin/*

doc:
	cd src && doxygen Doxyfile

$(OBJ_DIR)/jsoncpp/%.o: $(JSONCPP_SRC)/%.cpp
	test -d build/jsoncpp || mkdir -p build/jsoncpp
	$(CC) -c -o $@ $< $(CFLAGS)

$(OBJ_DIR)/%.o: src/%.cpp
	test -d build || mkdir build
	$(CC) -c -o $@ $< $(CFLAGS)

$(OBJ_DIR)/mcmc.o: src/mcmc.cpp
	test -d build || mkdir build
	$(CUDACC) -c -o $@ $< $(CFLAGS)

$(OBJ_DIR)/nll_kernels.o: src/nll_kernels.cpp
	test -d build || mkdir build
	$(CUDACC) -c -o $@ $< $(CFLAGS)

$(OBJ_DIR)/pdfz.o: src/pdfz.cpp
	test -d build || mkdir build
	$(CUDACC) -c -o $@ $< $(CFLAGS)

$(EXE): $(OBJECTS) $(JSONCPP_OBJECTS) $(OBJ_DIR)/mcmc.o $(OBJ_DIR)/nll_kernels.o $(OBJ_DIR)/pdfz.o
	test -d bin || mkdir bin
	$(GCC) -o $@ $^ $(CFLAGS) $(LFLAGS) $(CUDA_LFLAGS)


###### Test Infrastructure ############
test: bin/test_sxmc bin/bench_sxmc

GTEST_DIR = contrib/gtest-1.6.0
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# Generic Google Test code
build/gtest-all.o : $(GTEST_SRCS_)
	$(CXX) -I$(GTEST_DIR)/include -I$(GTEST_DIR) $(CXXFLAGS) -o $@ -c $(GTEST_DIR)/src/gtest-all.cc

build/gtest_main.o : $(GTEST_SRCS_)
	$(CXX) -I$(GTEST_DIR)/include -I$(GTEST_DIR) $(CXXFLAGS) -o $@ -c $(GTEST_DIR)/src/gtest_main.cc

build/gtest.a : build/gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

build/gtest_main.a : build/gtest-all.o build/gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

# Our test code
$(OBJ_DIR)/test/%.o: test/%.cpp
	test -d build/test || mkdir -p build/test
	$(CUDACC) -c -o $@ $< $(CFLAGS) -I$(GTEST_DIR)/include

bin/test_sxmc: $(TEST_OBJECTS) $(SXMC_NO_MAIN_FUNCTION_OBJECTS) build/gtest_main.a
	$(GCC) -o $@ $^ $(GCCFLAGS) $(LFLAGS) $(CUDA_LFLAGS)

## Benchmark code
$(OBJ_DIR)/bench_sxmc.o: bench/bench_sxmc.cpp
	$(CUDACC) -c -o $@ $< $(CFLAGS) -I$(GTEST_DIR)/include

bin/bench_sxmc: $(OBJ_DIR)/bench_sxmc.o $(SXMC_NO_MAIN_FUNCTION_OBJECTS)
	$(GCC) -o $@ $^ $(GCCFLAGS) $(LFLAGS) $(CUDA_LFLAGS)
