TARGETS = DAQStorage200s_Experiment

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

CXXFLAGS = -Wall -O2 $(ROOTFLAGS)
CXXLIBS = $(ROOTLIBS)

DAQStorageExperimet: DAQStorage200s_Experiment.o
	g++ -o $(TARGETS) DAQStorage200s_Experiment.o $(CXXLIBS) 
.cc.o:
	g++ -c $(CXXFLAGS) DAQStorage200s_Experiment.cc
