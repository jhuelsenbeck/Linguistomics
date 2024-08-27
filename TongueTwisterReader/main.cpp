#include <filesystem>
#include <iostream>
#include <string>
#include <fstream>
#include "McmcSummary.hpp"
#include "Msg.hpp"
#include "RandomVariable.hpp"
#include "UserSettings.hpp"

void printHeader(void);
void readConfigurationFile(std::string pathName, McmcSummary& summary);
void readDirectory(std::string filePath, int bi, McmcSummary& summary);
int readNumStates(std::string pathName, McmcSummary& summary);
Partition* readPartition(std::string pathName, McmcSummary& summary);



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
        }

    // read the third set of input, if available
    if (settings.getInputFile3() != "")
        {
        McmcSummary summary3(&rv);
        readConfigurationFile(settings.getInputFile3(), summary3);
        readDirectory(settings.getInputFile3(), settings.getBurnIn(), summary3);
        summary1 += summary3;
        }

    auto pathName = settings.getPath();
    auto findex = new std::ofstream(pathName + "/alignments.nytril", std::ios::out);
    summary1.output(pathName, *findex);
    findex->close();
    delete findex;
  
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

void readDirectory(std::string folder, int bi, McmcSummary& summary) {

    Partition* part = readPartition(folder, summary);
    if (part == NULL)
        Msg::warning("Could not find partition in configuration file");
        
    int numStates = readNumStates(folder, summary);
    std::cout << "  * Number of states = " << numStates << std::endl;
    if (numStates == 0)
        Msg::error("Could not find the NumberOfStates key in the JSON file");
    summary.setNumStates(numStates);

    for (const auto& entry : std::filesystem::directory_iterator(folder))
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
            summary.readAlnFile(filePath, bi, part);
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

int readNumStates(std::string pathName, McmcSummary& summary) {

    int ns = 0;
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
            {
            std::cout << "Reading number of states from " << filePath << std::endl;
            ns = summary.readNumStates(filePath);
            }
        }
    return ns;
}

Partition* readPartition(std::string pathName, McmcSummary& summary) {

    Partition* part = NULL;
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
            part = summary.readPartition(filePath);
        }
    return part;
}
