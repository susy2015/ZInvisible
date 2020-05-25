IDIR       = .
ODIR       = obj
SDIR       = .
RSDIR      = $(CMSSW_BASE)/src/SusyAnaTools/Tools

ZINVISIBLE = $(CMSSW_BASE)/src/ZInvisible

JSONDIR    = ${CMSSW_BASE}/src/json

# TODO for the brave: use Makefile.in, configure.ac, and configure script to create Makefile from Makefile.in
# Then you can use sytax like this:
# ZINVISIBLE = @ZINVISIBLE@ 

CXX        = g++

CXXFLAGS  += -I. -I$(CMSSW_BASE)/src -std=c++0x -I$(LHAPDF_DATA_PATH)/../../include
## Optimization flag
CXXFLAGS += -g #-O3
## Enable the maximun warning
#CXXFLAGS += -Wall -Wextra -Weffc++ -g

PYTHONCFG = python-config

## Include ROOT
CXXFLAGS  += $(shell root-config --cflags)

CXXDEPFLAGS = -MMD -MP

LD         = g++
LDFLAGS    =

INCLUDES += $(shell $(PYTHONCFG) --includes)

# include "include" directory
INCLUDES += -I$(ZINVISIBLE)/Tools/include -I$(JSONDIR)/include

CXXFLAGS += -DDOTENSORFLOWPYBIND

LIBS       = $(shell root-config --glibs) -lASImage
MT2LIB     = -L$(shell $(PYTHONCFG) --prefix)/lib $(shell $(PYTHONCFG) --libs)
LHAPDFLIB  = -L$(LHAPDF_DATA_PATH)/../../lib -lLHAPDF

PROGRAMS = makePlots

all: mkobj PyWraps $(PROGRAMS)


mkobj:
	@mkdir -p obj

PyWraps: $(ODIR)/samplesModule.so

$(ODIR)/samplesModule.so: $(ODIR)/samplesPyWrap.o $(ODIR)/samplesModulePyWrap.o
	$(CXX) -shared -o $@ $^

$(ODIR)/samplesPyWrap.o: $(RSDIR)/samples.cc $(RSDIR)/samples.h 
	$(CXX) --std=c++11 -c -fPIC -o $@ $<

$(ODIR)/samplesModulePyWrap.o: $(RSDIR)/samplesModule.cc
	$(CXX) --std=c++11 -c -fPIC -o $@ $<

$(ODIR)/%.o : $(SDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RSDIR) $(INCLUDES) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RSDIR) $(INCLUDES) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RSDIR) $(INCLUDES) -o $@ -c $<

$(ODIR)/%.o : $(RSDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RSDIR) $(INCLUDES) -o $@ -c $<

$(ODIR)/%.o : $(RSDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RSDIR) $(INCLUDES) -o $@ -c $<

$(ODIR)/%.o : $(RSDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RSDIR) $(INCLUDES) -o $@ -c $<

makePlots: $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/samples.o $(ODIR)/Plotter.o $(ODIR)/MakePlots.o $(ODIR)/PDFUncertainty.o $(ODIR)/RegisterFunctions.o $(ODIR)/SusyUtility.o $(ODIR)/customize.o
	$(LD) $^ $(LIBS) $(MT2LIB) $(LHAPDFLIB) -o $@

clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.so $(ODIR)/*.d $(PROGRAMS) core 

-include $(ODIR)/*.d
