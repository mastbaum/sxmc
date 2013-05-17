INCLUDE = -Isrc -I$(RATROOT)/include -I$(ROOTSYS)/include -I$(RATROOT)/src/stlplus -Icontrib/hemi -I/usr/local/cuda/include -I/opt/cuda-5.0/include
CFLAGS = -DVERBOSE=true -g $(INCLUDE) 
GCCFLAGS = -Wall -Werror -ffast-math -fdiagnostics-show-option
NVCCFLAGS = -arch=sm_30 -use_fast_math -lcudart
ROOTLIBS =  -lCore -lCint -lRIO -lMathCore -lHist -lGpad -lTree -lTree -lm
LFLAGS = -L$(RATROOT)/lib -lRATEvent_$(RATSYSTEM) -L$(ROOTSYS)/lib $(ROOTLIBS) -ljsoncpp

ifndef CUDA_ROOT
$(warning *** CUDA_ROOT is not set, defaulting to CPU-only build ***)
CC = g++ $(GCCFLAGS) -DHEMI_CUDA_DISABLE
CUDACC = $(CC)
else
CC = $(CUDA_ROOT)/bin/nvcc $(NVCCFLAGS) -I$(CUDA_ROOT)/include
CUDACC = $(CC) -x cu
CC += --compiler-options "$(GCCFLAGS) -Wno-unused-function"
endif

OBJ_DIR = build
SOURCES = $(filter-out src/nll.cpp, $(wildcard src/*.cpp))
OBJECTS = $(SOURCES:src/%.cpp=$(OBJ_DIR)/%.o)
EXE = sensitivity

ifndef RATROOT
$(error RATROOT is not set)
endif

ifndef ROOTSYS
$(error ROOTSYS is not set)
endif

all: $(OBJ_DIR)/nll.o $(OBJECTS) $(EXE)

.PHONY: doc

clean:
	-$(RM) build/*.o bin/*

doc:
	cd src && doxygen Doxyfile

$(OBJ_DIR)/%.o: src/%.cpp
	test -d build || mkdir build
	$(CC) -c -o $@ $< $(CFLAGS) $(LFLAGS)

$(OBJ_DIR)/nll.o: src/nll.cpp
	test -d build || mkdir build
	$(CUDACC) -c -o $@ $< $(CFLAGS) $(LFLAGS)

$(EXE): $(OBJECTS) $(CUDA_OBJECTS) $(OBJ_DIR)/nll.o
	test -d bin || mkdir bin
	$(CC) -o bin/$@ $^ $(CFLAGS) $(LFLAGS) $(CUDA_LFLAGS)

