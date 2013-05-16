CUDA_INCLUDE = -I$(CUDA_ROOT)/include
CUDA_CFLAGS = -arch=sm_30 -O3 -use_fast_math -I$(ROOTSYS)/include
CUDA_LFLAGS = -L$(CUDA_ROOT)/lib64 -lcudart -L$(ROOTSYS)/lib -lCore -lCint -lRIO -lMathCore -lm
INCLUDE = -Isrc -I$(RATROOT)/include -I$(ROOTSYS)/include -I$(RATROOT)/src/stlplus
CFLAGS = -DVERBOSE=true -O3 $(INCLUDE) 
GCCFLAGS = --std=c++0x -Wall -Werror -ffast-math -fdiagnostics-show-option
LFLAGS = -L$(RATROOT)/lib -lRATEvent_$(RATSYSTEM) $(shell root-config --libs) -ljsoncpp

CC = g++
NVCC = $(CUDA_ROOT)/bin/nvcc

OBJ_DIR = build
SOURCES = $(wildcard src/*.cpp)
CUDA_SOURCES = $(wildcard src/*.cu)
OBJECTS = $(SOURCES:src/%.cpp=$(OBJ_DIR)/%.o)
CUDA_OBJECTS = $(CUDA_SOURCES:src/%.cu=$(OBJ_DIR)/%.o)
EXE = sensitivity

ifndef RATROOT
$(error RATROOT is not set)
endif

ifndef ROOTSYS
$(error ROOTSYS is not set)
endif

ifndef CUDA_ROOT
$(error CUDA_ROOT is not set)
endif

all: $(CUDA_OBJECTS) $(OBJECTS) $(EXE)

clean:
	-$(RM) build/*.o bin/*

$(OBJ_DIR)/%.o: src/%.cpp
	test -d build || mkdir build
	$(CC) -c -o $@ $< $(GCCFLAGS) $(CFLAGS) $(LFLAGS)

$(OBJ_DIR)/nll.o:
	test -d build || mkdir build
	$(NVCC) -c -o build/nll.o src/nll.cu $(CUDA_CFLAGS) $(CUDA_LFLAGS)

$(EXE): $(OBJECTS) $(CUDA_OBJECTS)
	test -d bin || mkdir bin
	$(CC) -o bin/$@ $^ $(CFLAGS) $(LFLAGS) $(CUDA_LFLAGS)

