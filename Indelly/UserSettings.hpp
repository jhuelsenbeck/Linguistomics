#ifndef UserSettings_hpp
#define UserSettings_hpp

#include <string>
#include "json.hpp"


enum SubstitutionModel { jc69, gtr, custom };

class UserSettings {

    public:
        static UserSettings&    userSettings(void) {
                                    static UserSettings us;
                                    return us;
                                }
        bool                    getCalculateMarginalLikelihood(void) { return calculateMarginalLikelihood; }
        int                     getCheckPointFrequency(void) { return checkPointFrequency; }
        std::string             getDataFile(void) { return dataFile; }
        double                  getBranchLengthLambda(void) { return branchLengthLambda; }
        int                     getBurninLength(void) { return burninLength; }
        std::string             getOutFile(void) { return outFile; }
        int                     getNumChains(void) { return numChains; }
        int                     getNumMcmcCycles(void) { return numMcmcCycles; }
        int                     getNumIndelCategories(void) { return numIndelCategories; }
        int                     getNumRateCategories(void) { return numRateCategories; }
        int                     getNumTunes(void) { return numTunes; }
        int                     getPreburninLength(void) { return preburninLength; }
        int                     getPrintFrequency(void) { return printFrequency; }
        int                     getSampleFrequency(void) { return sampleFrequency; }
        int                     getSampleLength(void) { return sampleLength; }
        int                     getSampleToStoneFrequency(void) { return sampleToStoneFrequency; }
        uint32_t                getSeed(void) { return seed; }
        int                     getSubstitutionModel(void) { return substitutionModel; }
        double                  getTemperature(void) { return temperature; }
        int                     getTuneLength(void) { return tuneLength; }
        void                    print(void);
        void                    readCommandLineArguments(int argc, char* argv[]);
        bool                    getUseOnlyCompleteWords(void) { return useOnlyCompleteWords; }
    
    private:
                                UserSettings(void);
                                UserSettings(const UserSettings& s) = delete;
        void                    usage(void);
        std::string             getVariable(nlohmann::json& settings, const char *name);
        std::string             outFile;
        std::string             dataFile;
        std::string             executablePath;
        double                  branchLengthLambda;
        int                     checkPointFrequency;
        int                     numIndelCategories;
        int                     numMcmcCycles;
        int                     numRateCategories;
        int                     printFrequency;
        int                     sampleFrequency;
        int                     substitutionModel;
        bool                    useOnlyCompleteWords;
        bool                    calculateMarginalLikelihood;
        int                     preburninLength;
        int                     numTunes;
        int                     tuneLength;
        int                     burninLength;
        int                     sampleLength;
        int                     sampleToStoneFrequency;
        uint32_t                seed;
        int                     numChains;
        double                  temperature;
};

#endif
