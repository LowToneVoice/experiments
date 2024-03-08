TARGETS = Various_ImportTime

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

CXXFLAGS = -Wall -O2 $(ROOTFLAGS)
CXXLIBS = $(ROOTLIBS)

DAQStorageExperimet: Various_ImportTime.o
	g++ -o $(TARGETS) Various_ImportTime.o $(CXXLIBS) 
.cc.o:
	g++ -c $(CXXFLAGS) Various_ImportTime.cc
