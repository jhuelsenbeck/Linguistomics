#include <filesystem>
#include <iostream>
#include <string>
#include "McmcSummary.hpp"
#include "RandomVariable.hpp"
#include "UserSettings.hpp"

void printHeader(void);
void readConfigurationFile(UserSettings& settings, McmcSummary& summary);



int main(int argc, char* argv[]) {

    // TODO: Output information on partition group frequencies (i.e., size of circles), with CIs.
    //       Check that the partition group indices are consistent across all output in this program (esp. check .csv file)

    RandomVariable rv;
    printHeader();

    // get the user settings
    UserSettings& settings = UserSettings::userSettings();
    settings.readCommandLineArguments(argc, argv);
    settings.print();
    
    // go through all of the files in the directory
    McmcSummary summary(&rv);
    readConfigurationFile(settings, summary); // this will instantiate the partition sets object, if partitions are present
    for (const auto& entry : std::filesystem::directory_iterator(settings.getPath()))
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
            summary.readTsvFile(filePath, settings.getBurnIn());
            }
        else if (fileExtension == ".aln")
            {
            show = true;
            summary.readAlnFile(filePath, settings.getBurnIn());
            }
        else if (fileExtension == ".tre")
            {
            show = true;
            summary.readTreFile(filePath, settings.getBurnIn());
            }

        show = false;
        if (show)
            std::cout << filePath << std::endl;
        }
    
    //summary.print();
    summary.output(settings);

    return 0;
}

void printHeader(void) {

    std::cout << "   TongueTwister Reader" << std::endl;
    std::cout << "   * John P. Huelsenbeck" << std::endl;
}

void readConfigurationFile(UserSettings& settings, McmcSummary& summary) {

    for (const auto& entry : std::filesystem::directory_iterator(settings.getPath()))
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
            summary.readConfigFile(fp);
        }
}
