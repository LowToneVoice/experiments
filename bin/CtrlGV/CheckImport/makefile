TARGETS = CheckImportTime

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

CXXFLAGS = -Wall -O2 $(ROOTFLAGS)
CXXLIBS = $(ROOTLIBS)

DAQStorageExperimet: CheckImportTime.o
	g++ -o $(TARGETS) CheckImportTime.o $(CXXLIBS) 
.cc.o:
	g++ -c $(CXXFLAGS) CheckImportTime.cc
