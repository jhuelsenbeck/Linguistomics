#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include "Msg.hpp"
#include "UserSettings.hpp"



UserSettings::UserSettings(void) {

    // dafault values
    inFile  = "";
    outFile = "";
    burnIn  = 0;
}

void UserSettings::readCommandLineArguments(int argc, char* argv[]) {

    // move the command-line arguments to a good-ol' C++ STL vector of strings
    std::vector<std::string> commands;
    for (int i=0; i<argc; i++)
        commands.push_back(argv[i]);
    
    // check that we have an odd number of commands
    if (commands.size() % 2 == 0 || commands.size() == 1)
        {
        usage();
        Msg::error("Problem with arguments");
        }
    
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
            if (arg == "-i")
                inFile = cmd;
            else if (arg == "-o")
                outFile = cmd;
            else if (arg == "-b")
                burnIn = atoi(cmd.c_str());
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
    std::cout << "   * Input path name  = \"" << inFile << "\"" << std::endl;
    std::cout << "   * Output file name = \"" << outFile << "\"" << std::endl;
    std::cout << "   * Burn in number   = " << burnIn << std::endl;
    std::cout << std::endl;
}

void UserSettings::usage(void) {

    std::cout << "   Usage" << std::endl;
    std::cout << "   * -i -- File with data to be summarized" << std::endl;
    std::cout << "   * -o -- File name for output" << std::endl;
    std::cout << "   * -b -- Burn-in number" << std::endl;
    std::cout << std::endl;
}

