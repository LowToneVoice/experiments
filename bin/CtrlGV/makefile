#target
#TARGET = DAQStorageExperiment_2021
TARGET = ValveControll_EDMcellExperiment
#TARGET = ValveControll_EDMcellExperimentValveClosed
#TARGET = ValveControll_EDMFillSearch

HEADER = settings

#suffix
ObjSuf          = .o
SrcSuf          = .cc
ExeSuf          = 
HdrSuf          = .h

#set suffix
TARGETOBJ = $(TARGET)$(ObjSuf)
TARGETSRC = $(TARGET)$(SrcSuf)
TARGETEXE = $(TARGET)$(ExeSuf)
TARGETHDR = $(HEADER)$(HdrSuf)

#for ROOT
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)

#for linux
CXX       = g++
CXXFLAGS  = -O2 -Wall
LD        = g++
LDFLAGS   = -O2 -Wall

# create parameters
CXXFLAGS += $(ROOTCFLAGS)
LIBS      = $(ROOTLIBS) $(SYSLIBS)
GLIBS     = $(ROOTGLIBS) $(SYSLIBS)

all: $(TARGETEXE)

$(TARGETEXE): $(TARGETOBJ)
	$(LD) $(LDFLAGS) $(TARGETOBJ) $(LIBS) -o $(TARGETEXE)

$(TARGETOBJ): $(TARGETSRC)
	$(CXX) $(CXXFLAGS) -c -o $(TARGETOBJ) $(TARGETSRC)

$(TARGETOBJ): $(TARGETHDR)

clean:
	rm -f *.o
	rm -f *~
