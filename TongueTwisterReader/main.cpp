#include <filesystem>
#include <iostream>
#include <string>
#include "McmcSummary.hpp"
#include "RandomVariable.hpp"
#include "UserSettings.hpp"

void printHeader(void);
void readConfigurationFile(std::string pathName, McmcSummary& summary);
void readDirectory(std::string filePath, int bi, McmcSummary& summary);



int main(int argc, char* argv[]) {

    RandomVariable rv;
    printHeader();

    // get the user settings
    UserSettings& settings = UserSettings::userSettings();
    settings.readCommandLineArguments(argc, argv);
    settings.print();
    
    // go through all of the files in the directory
    McmcSummary summary1(&rv);
    readConfigurationFile(settings.getInputFile(), summary1); // this will instantiate the partition sets object, if partitions are present
    readDirectory(settings.getInputFile(), settings.getBurnIn(), summary1);
    
    // read the second set of input, if available
    if (settings.getInputFile2() != "")
        {
        McmcSummary summary2(&rv);
        readConfigurationFile(settings.getInputFile2(), summary2);
        readDirectory(settings.getInputFile2(), settings.getBurnIn(), summary2);
        summary1 += summary2;
        summary1.output(settings.getPath());
        }
    else
        {
        summary1.output(settings.getPath());
        }

    return 0;
}

void printHeader(void) {

    std::cout << "   TongueTwister Reader" << std::endl;
    std::cout << "   * John P. Huelsenbeck" << std::endl;
}

void readConfigurationFile(std::string pathName, McmcSummary& summary) {

    for (const auto& entry : std::filesystem::directory_iterator(pathName))
        {
        std::filesystem::path fp = entry.path();
        #ifdef _CONSOLE
        auto filePath      = fp.string();
        auto fileExtension = fp.extension();
        #else
        auto filePath      = fp;
        auto fileExtension = fp.extension();
        #endif

        if (fileExtension == ".config")
            summary.readConfigFile(filePath);
        }
}

void readDirectory(std::string filePath, int bi, McmcSummary& summary) {

    for (const auto& entry : std::filesystem::directory_iterator(filePath))
        {
        std::filesystem::path fp = entry.path();
        #ifdef _CONSOLE
        auto filePath      = fp.string();
        auto fileExtension = fp.extension();
        #else
        auto filePath      = fp;
        auto fileExtension = fp.extension();
        #endif
        bool show = false;
                
        if (fileExtension == ".tsv")
            {
            show = true;
            summary.readTsvFile(filePath, bi);
            }
        else if (fileExtension == ".aln")
            {
            show = true;
            summary.readAlnFile(filePath, bi);
            }
        else if (fileExtension == ".tre")
            {
            show = true;
            summary.readTreFile(filePath, bi);
            }

        if (show)
            std::cout << filePath << std::endl;
        }
}
