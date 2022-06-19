CC = g++
CXX = g++

CFLAGS = -std=c++20 -O3 -fopenmp $(INCLUDES)
CXXFLAGS = -std=c++20 -O3 -fopenmp $(INCLUDES)

INCLUDES =
LDFLAGS =
LDLIBS = 

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
