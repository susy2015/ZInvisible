IDIR       = .
ODIR       = obj
SDIR       = .
RSDIR      = $(CMSSW_BASE)/src/SusyAnaTools/Tools

CXX        = g++

CXXFLAGS  += -I. -I$(CMSSW_BASE)/src -std=c++0x
## Optimization flag
CXXFLAGS += -g #-O3
## Enable the maximun warning
#CXXFLAGS += -Wall -Wextra -Weffc++ -g

## Include ROOT
CXXFLAGS  += $(shell root-config --cflags)

CXXDEPFLAGS = -MMD -MP

LD         = g++
LDFLAGS    =

LIBS       = $(shell root-config --glibs)
MT2LIB     = -L$(CMSSW_BASE)/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2

#OBJS       = $(patsubst %, $(ODIR)/%, $(OBJ))
PROGRAMS = calcEff harvestHists makePlots nEvts

all: $(PROGRAMS)


$(ODIR)/%.o : $(SDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(RSDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(RSDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

calcEff: $(ODIR)/NTupleReader.o $(ODIR)/samples.o $(ODIR)/calcEff.o
	$(LD) $^ $(LIBS) $(MT2LIB) -o $@

harvestHists: $(ODIR)/harvestHists.o
	$(LD) $^ $(LIBS) -o $@

makePlots: $(ODIR)/NTupleReader.o $(ODIR)/samples.o $(ODIR)/Plotter.o $(ODIR)/MakePlots.o
	$(LD) $^ $(LIBS) $(MT2LIB) -o $@

nEvts: $(ODIR)/samples.o $(ODIR)/nEvts.o
	$(LD) $^ $(LIBS) -o $@

clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.d $(PROGRAMS) core 

-include $(ODIR)/*.d
