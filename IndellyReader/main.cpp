#include <filesystem>
#include <iostream>
#include <string>
#include "McmcSummary.hpp"
#include "RandomVariable.hpp"
#include "UserSettings.hpp"

void printHeader(void);



int main(int argc, char* argv[]) {

    RandomVariable rv;
    printHeader();

    // get the user settings
    UserSettings& settings = UserSettings::userSettings();
    settings.readCommandLineArguments(argc, argv);
    settings.print();
    
    // go through all of the files in the directory
    McmcSummary summary(&rv);
    for (const auto& entry : std::filesystem::directory_iterator(settings.getPath()))
        {
        std::filesystem::path filePath = entry.path();
        std::filesystem::path fileName = filePath.filename();
        std::filesystem::path fileExtension = filePath.extension();
        std::cout << filePath << std::endl;
        
        if (fileExtension == ".tsv")
            {
            summary.readTsvFile(filePath, settings.getBurnIn());
            }
        else if (fileExtension == ".aln")
            {
            summary.readAlnFile(filePath, settings.getBurnIn());
            }
        else if (fileExtension == ".tre")
            {
            summary.readTreFile(filePath, settings.getBurnIn());
            }
        }
    
    summary.print();


    return 0;
}

void printHeader(void) {

    std::cout << "   TongueTwister Reader" << std::endl;
    std::cout << "   * John P. Huelsenbeck" << std::endl;
}
