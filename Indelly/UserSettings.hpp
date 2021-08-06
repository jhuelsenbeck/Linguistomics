#ifndef UserSettings_hpp
#define UserSettings_hpp

#include <string>



class UserSettings {

    public:
        static UserSettings&    userSettings(void) {
                                    static UserSettings us;
                                    return us;
                                }
        std::string             getDataFile(void) { return dataFile; }
        double                  getInverseTreeLength(void) { return inverseTreeLength; }
        std::string             getOutFile(void) { return outFile; }
        int                     getNumMcmcCycles(void) { return numMcmcCycles; }
        int                     getPrintFrequency(void) { return printFrequency; }
        int                     getSampleFrequency(void) { return sampleFrequency; }
        std::string             getSubstitutionModel(void) { return substitutionModel; }
        void                    print(void);
        void                    readCommandLineArguments(int argc, char* argv[]);
    
    private:
                                UserSettings(void);
                                UserSettings(const UserSettings& s) = delete;
        void                    usage(void);
        std::string             dataFile;
        std::string             outFile;
        int                     numMcmcCycles;
        int                     printFrequency;
        int                     sampleFrequency;
        double                  inverseTreeLength;
        std::string             substitutionModel;
};

#endif
