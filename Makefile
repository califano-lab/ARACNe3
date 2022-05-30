CC = g++
CXX = g++

CXXFLAGS = -std=c++20 $(INCLUDES)

INCLUDES = -I/opt/homebrew/lib/R/4.1/site-library/Rcpp/include/ -I/Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/include/
LDFLAGS = -L/opt/homebrew/Cellar/boost/1.78.0_1/lib/
LDLIBS = 


ARACNe3: APMI.o MatrixReglistIO.o

NullModel: APMI.o

APMI: 

MatrixReglistIO: 

APMI.o:

MatrixReglistIO.o:

.PHONY: all
all: ARACNe3 APMI MatrixReglistIO

.PHONY: clean
clean: 
	rm *.o APMI MatrixReglistIO NullModel
