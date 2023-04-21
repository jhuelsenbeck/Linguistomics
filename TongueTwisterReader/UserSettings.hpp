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
        std::string             getInputFile(void) { return inFile; }
        std::string             getInputFile2(void) { return inFile2; }
        std::string             getInputFile3(void) { return inFile3; }
        std::string             getOutFile(void) { return outFile; }
        int                     getBurnIn(void) { return burnIn; }
        std::string             getPath(void) { return inFile; }
        void                    print(void);
        void                    readCommandLineArguments(int argc, char* argv[]);
    
    private:
                                UserSettings(void);
                                UserSettings(const UserSettings& s) = delete;
        void                    usage(void);
        std::string             inFile;
        std::string             inFile2;
        std::string             inFile3;
        std::string             outFile;
        int                     burnIn;
};

#endif
