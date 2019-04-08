CXX=g++
CXXFLAGS=-Wall -O3

TESTFLAGS=-DNSLICES=10 -DNRADIANS=11 

# Uncomment the following to enable report of timings
#EXTRA+=-D_PRINT_MYTIMER_

HEADER_FILES=DefsMTUsed.h Elem_Stiff_Mat.h FibrilIdentifier.h InputData.h Loading.h MeshingData.h PostProcess.h SolveSystem.h MyTimer.h

all: test production

test: ./tendonmech_test
production: ./tendonmech

./tendonmech: main.cpp $(HEADER_FILES)
	$(CXX) $(CXXFLAGS) -o ./tendonmech main.cpp

./tendonmech_test: main.cpp $(HEADER_FILES)
	$(CXX) $(CXXFLAGS) $(TESTFLAGS) -o ./tendonmech_test main.cpp

clean:
	rm -f ./tendonmech ./tendonmech_test

.PHONY: test production

