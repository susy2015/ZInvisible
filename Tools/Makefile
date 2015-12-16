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
PROGRAMS = calcEff harvestHists makePlots systematics nJetWgtSyst

all: mkobj sampPyWrap $(PROGRAMS)


mkobj:
	@mkdir -p obj

sampPyWrap: $(ODIR)/samplesModule.so

$(ODIR)/samplesModule.so: $(ODIR)/samplesPyWrap.o $(ODIR)/samplesModulePyWrap.o
	$(CXX) -shared -o $@ $^

$(ODIR)/samplesPyWrap.o: $(RSDIR)/samples.cc $(RSDIR)/samples.h 
	$(CXX) --std=c++11 -c -fPIC -o $@ $<

$(ODIR)/samplesModulePyWrap.o: $(RSDIR)/samplesModule.cc
	$(CXX) --std=c++11 -c -fPIC -o $@ $<


$(ODIR)/%.o : $(SDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RSDIR) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RSDIR) -o $@ -c $<

$(ODIR)/%.o : $(RSDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RSDIR) -o $@ -c $<

$(ODIR)/%.o : $(RSDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RSDIR) -o $@ -c $<

calcEff: $(ODIR)/NTupleReader.o $(ODIR)/samples.o $(ODIR)/calcEff.o $(ODIR)/RegisterFunctions.o $(ODIR)/Systematic.o $(ODIR)/Plotter.o $(ODIR)/MiniTupleMaker.o
	$(LD) $^ $(LIBS) $(MT2LIB) -o $@

harvestHists: $(ODIR)/harvestHists.o
	$(LD) $^ $(LIBS) -o $@

makePlots: $(ODIR)/NTupleReader.o $(ODIR)/samples.o $(ODIR)/Plotter.o $(ODIR)/MakePlots.o $(ODIR)/MiniTupleMaker.o $(ODIR)/RegisterFunctions.o $(ODIR)/Systematic.o 
	$(LD) $^ $(LIBS) $(MT2LIB) -o $@

systematics: $(ODIR)/NTupleReader.o $(ODIR)/samples.o $(ODIR)/Systematic.o  $(ODIR)/Plotter.o $(ODIR)/systematics.o $(ODIR)/MiniTupleMaker.o $(ODIR)/RegisterFunctions.o
	$(LD) $^ $(LIBS) $(MT2LIB) -o $@

nJetWgtSyst: $(ODIR)/NTupleReader.o $(ODIR)/nJetWgtSyst.o $(ODIR)/samples.o
	$(LD) $^ $(LIBS) $(MT2LIB) -o $@

clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.so $(ODIR)/*.d $(PROGRAMS) core 

-include $(ODIR)/*.d
