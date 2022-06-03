CC = g++
CXX = g++

CXXFLAGS = -std=c++20 $(INCLUDES)

INCLUDES =
LDFLAGS =
LDLIBS = 

.PHONY: release
release: ARACNe3 clean

ARACNe3: NullModel.o MatrixReglistIO.o APMI.o FDRPruning.o MaxEntPruning.o RegWebFns.o

NullModel.o:

MatrixReglistIO.o:

APMI.o:

FDRPruning.o:

MaxEntPruning.o:

RegWebFns.o:

.PHONY: all
all: ARACNe3 APMI MatrixReglistIO

.PHONY: clean
clean: 
	rm *.o APMI MatrixReglistIO NullModel
