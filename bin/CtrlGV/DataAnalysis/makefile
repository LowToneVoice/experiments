TARGETS = MakeRootForUCNStorage

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

CXXFLAGS = -Wall -O2 $(ROOTFLAGS)
CXXLIBS = $(ROOTLIBS)

DAQStorageExperimet: MakeRootForUCNStorage.o
	g++ -o $(TARGETS) MakeRootForUCNStorage.o $(CXXLIBS) 
.cc.o:
	g++ -c $(CXXFLAGS) MakeRootForUCNStorage.cc
