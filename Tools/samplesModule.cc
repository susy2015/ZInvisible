#include "../../SusyAnaTools/Tools/samples.h"

#include <string>

extern "C" {
    AnaSamples::SampleCollection* SC_new()
    {
        AnaSamples::SampleSet *ss = new AnaSamples::SampleSet();
        return new AnaSamples::SampleCollection(*ss); 
    }
    int SC_samples_size(AnaSamples::SampleCollection* sc, char *scn){ return (*sc)[std::string(scn)].size(); }
    char const ** SC_samples(AnaSamples::SampleCollection* sc, char *scn)
    {
        auto& sampleVec = (*sc)[std::string(scn)];
        const char **array = new const char*[sampleVec.size()];
        int i = 0;
        for(auto& sample : sampleVec)
        {
            array[i++] = sample.filePath.c_str();
        }
        return array;
    }
    char const ** SC_samples_names(AnaSamples::SampleCollection* sc, char *scn)
    {
        auto& sampleVec = sc->getSampleLabels(std::string(scn));
        const char **array = new const char*[sampleVec.size()];
        int i = 0;
        for(auto& sample : sampleVec)
        {
            array[i++] = sample.c_str();
        }
        return array;
    }
    int SC_samplecollection_size(AnaSamples::SampleCollection* sc, char *scn){ return sc->size(); }
    char const ** SC_samplecollection_names(AnaSamples::SampleCollection* sc)
    {
        const char **array = new const char*[sc->size()];
        int i = 0;
        for(auto& sample : *sc)
        {
            array[i++] = sample.first.c_str();
        }
        return array;
    }
}
