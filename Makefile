CXX=g++
CXXFLAGS=-g -Wall -std=c++11
PROG=TestQR
SRC=QR.cpp Example.cpp

all: $(PROG)

$(PROG) : $(SRC)
	$(CXX) $(CXXFLAGS) -o $(PROG) $(SRC)

clean : $(PROG)
	rm $(PROG)