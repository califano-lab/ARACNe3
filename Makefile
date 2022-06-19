CC = clang++
CXX = clang++

CFLAGS = -std=c++20 -Xpreprocessor -fopenmp -O3 $(INCLUDES)
CXXFLAGS = -std=c++20 -Xpreprocessor -fopenmp -O3 $(INCLUDES)

INCLUDES = 
LDFLAGS = -L/opt/homebrew/opt/llvm/lib
LDLIBS = -lomp

.PHONY: release
release: ARACNe3 clean

ARACNe3: NullModel.o MatrixReglistIO.o APMI.o AlphaPruning.o MaxEntPruning.o RegWebFns.o Consolidator.o

NullModel.o:

MatrixReglistIO.o:

APMI.o:

AlphaPruning.o:

MaxEntPruning.o:

RegWebFns.o:

Consolidator.o:

.PHONY: all
all: ARACNe3 

.PHONY: clean
clean: 
	rm *.o *.dSYM
