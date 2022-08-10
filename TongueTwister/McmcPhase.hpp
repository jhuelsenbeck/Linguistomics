#ifndef McmcPhase_hpp
#define McmcPhase_hpp

#include <string>
#include <vector>


class McmcPhase {

    public:
                                    McmcPhase(void) = delete;
                                    McmcPhase(int pbl, int tl, int bl, int sl, int sf);
        int                         getPreBurnLength(void) { return preBurnLength; }
        int                         getBurnLength(void) { return burnLength; }
        int                         getTuneLength(void) { return tuneLength; }
        int                         getSampleLength(void) { return sampleLength; }
        int                         getSampleFrequency(void) { return sampleFrequency; }
        int                         getLengthForPhase(std::string ph);
        std::vector<std::string>&   getPhases(void) { return phases; }

    private:
        int                         preBurnLength;
        int                         tuneLength;
        int                         burnLength;
        int                         sampleLength;
        int                         sampleFrequency;
        std::vector<std::string>    phases;
};

#endif
