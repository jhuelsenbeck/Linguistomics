#include <iostream>
#include <vector>
#include "Msg.hpp"
#include "UserSettings.hpp"

#define DEBUG_MODE



UserSettings::UserSettings(void) {

    // dafault values
    dataFile          = "";
    treeFile          = "";
    numMcmcCycles     = 1000;
    printFrequency    = 100;
    sampleFrequency   = 100;
    inverseTreeLength = 1.0;
}

void UserSettings::readCommandLineArguments(int argc, char* argv[]) {

    // move the command-line arguments to a good-ol' C++ STL vector of strings
    std::vector<std::string> commands;
#   if defined(DEBUG_MODE)
    commands.push_back(argv[0]);
    commands.push_back("-d");
    commands.push_back("/Users/johnh/Desktop/Indelly/alignments_old");
    commands.push_back("-t");
    commands.push_back("/Users/johnh/Desktop/Indelly/words9.tree");
    commands.push_back("-o");
    commands.push_back("/Users/johnh/Desktop/Indelly/out/test");
    commands.push_back("-n");
    commands.push_back("100000");
    commands.push_back("-p");
    commands.push_back("1");
    commands.push_back("-s");
    commands.push_back("500");
    commands.push_back("-l");
    commands.push_back("0.15");
#   else
    for (int i=0; i<argc; i++)
        commands.push_back(argv[i]);
#   endif
    
    // read the command-line arguments, starting with the second one
    std::string arg = "";
    for (int i=1; i<commands.size(); i++)
        {
        std::string cmd = commands[i];
        if (arg == "")
            {
            arg = cmd;
            }
        else
            {
            if (arg == "-d")
                dataFile = cmd;
            else if (arg == "-t")
                treeFile = cmd;
            else if (arg == "-o")
                outFile = cmd;
            else if (arg == "-n")
                numMcmcCycles = atoi(cmd.c_str());
            else if (arg == "-p")
                printFrequency = atoi(cmd.c_str());
            else if (arg == "-s")
                sampleFrequency = atoi(cmd.c_str());
            else if (arg == "-l")
                inverseTreeLength = atof(cmd.c_str());
            else
                {
                usage();
                Msg::error("Unknown command \"" + arg + "\"");
                }
            
            arg = "";
            }
        }
}

void UserSettings::print(void) {

    std::cout << "   Settings" << std::endl;
    std::cout << "   * Directory with initial word alignments = \"" << dataFile << "\"" << std::endl;
    std::cout << "   * Tree file path and name                = \"" << treeFile << "\"" << std::endl;
    std::cout << "   * Number of MCMC cycles                  = " << numMcmcCycles << std::endl;
    std::cout << "   * Print-to-screen frequency              = " << printFrequency << std::endl;
    std::cout << "   * Chain sample frequency                 = " << sampleFrequency << std::endl;
    std::cout << "   * Tree length prior parameter            = " << inverseTreeLength << std::endl;
    std::cout << std::endl;
}

void UserSettings::usage(void) {

    std::cout << "   Usage" << std::endl;
    std::cout << "   * -d -- Path to initial alignments of words" << std::endl;
    std::cout << "   * -t -- Tree file path and name" << std::endl;
    std::cout << "   * -n -- Number of MCMC cycles" << std::endl;
    std::cout << "   * -p -- Print-to-screen frequency" << std::endl;
    std::cout << "   * -s -- Chain sample frequency" << std::endl;
    std::cout << "   * -l -- Inverse of the tree length prior" << std::endl;
    std::cout << std::endl;
}
