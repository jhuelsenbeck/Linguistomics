#ifndef UserSettings_hpp
#define UserSettings_hpp

#include <string>


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
        double                  getInverseTreeLength(void) { return inverseTreeLength; }
        std::string             getOutFile(void) { return outFile; }
        int                     getNumMcmcCycles(void) { return numMcmcCycles; }
        int                     getNumIndelCategories(void) { return numIndelCategories; }
        int                     getNumRateCategories(void) { return numRateCategories; }
        int                     getPrintFrequency(void) { return printFrequency; }
        int                     getSampleFrequency(void) { return sampleFrequency; }
        int                     getSubstitutionModel(void) { return substitutionModel; }
        bool                    getUseEigenSystem(void) { return useEigenSystem; }
        void                    print(void);
        void                    readCommandLineArguments(int argc, char* argv[]);
        bool                    getUseOnlyCompleteWords(void) { return useOnlyCompleteWords; }
    
    private:
                                UserSettings(void);
                                UserSettings(const UserSettings& s) = delete;
        void                    usage(void);
        bool                    calculateMarginalLikelihood;
        int                     checkPointFrequency;
        std::string             dataFile;
        std::string             executablePath;
        double                  inverseTreeLength;
        std::string             outFile;
        int                     numIndelCategories;
        int                     numMcmcCycles;
        int                     numRateCategories;
        int                     printFrequency;
        int                     sampleFrequency;
        int                     substitutionModel;
        bool                    useEigenSystem;
        bool                    useOnlyCompleteWords;
};

#endif
