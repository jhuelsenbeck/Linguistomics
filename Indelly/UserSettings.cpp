#include <iostream>
#include <vector>
#include "Msg.hpp"
#include "UserSettings.hpp"

#define DEBUG_MODE



UserSettings::UserSettings(void) {

    // dafault values
    dataFile                    = "";
    outFile                     = "";
    numMcmcCycles               = 1000;
    printFrequency              = 100;
    sampleFrequency             = 100;
    inverseTreeLength           = 1.0;
    substitutionModel           = jc69;
    calculateMarginalLikelihood = false;
}

void UserSettings::readCommandLineArguments(int argc, char* argv[]) {

    // move the command-line arguments to a good-ol' C++ STL vector of strings
    std::vector<std::string> commands;
    executablePath = argv[0];
#   if defined(DEBUG_MODE)
    commands.push_back(executablePath);
    commands.push_back("-d");
    commands.push_back("/Users/johnh/Repositories/Linguistomics/Run/Run1/config.json");
    commands.push_back("-o");
    commands.push_back("/Users/johnh/Desktop/Indelly/out/test_jc");
    commands.push_back("-m");
    commands.push_back("jc69");
    commands.push_back("-n");
    commands.push_back("100000");
    commands.push_back("-p");
    commands.push_back("100");
    commands.push_back("-s");
    commands.push_back("100");
    commands.push_back("-l");
    commands.push_back("0.15");
    commands.push_back("-z");
    commands.push_back("no");
#   else
    for (int i=0; i<argc; i++)
        commands.push_back(argv[i]);
#   endif

    if (commands.size() % 2 == 0 || commands.size() == 1)
        {
        usage();
        Msg::error("Problem with arguments for \"" + executablePath + "\"");
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
            if (arg == "-d")
                dataFile = cmd;
            else if (arg == "-o")
                outFile = cmd;
            else if (arg == "-m")
                {
                if (cmd == "jc69")
                    substitutionModel = jc69;
                else if (cmd == "gtr")
                    substitutionModel = gtr;
                else if (cmd == "custom")
                    substitutionModel = custom;
                else
                    Msg::error("Unknon substitution model " + cmd);
                }
            else if (arg == "-n")
                numMcmcCycles = atoi(cmd.c_str());
            else if (arg == "-p")
                printFrequency = atoi(cmd.c_str());
            else if (arg == "-s")
                sampleFrequency = atoi(cmd.c_str());
            else if (arg == "-l")
                inverseTreeLength = atof(cmd.c_str());
            else if (arg == "-z")
                {
                if (cmd == "yes")
                    calculateMarginalLikelihood = true;
                else if (cmd == "no")
                    calculateMarginalLikelihood = false;
                else
                    Msg::error("Unknon option for calculating marginal likelihood");
                }
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
    std::cout << "   * Executable path                        = " << executablePath << std::endl;
    std::cout << "   * File with initial word alignments      = \"" << dataFile << "\"" << std::endl;
    if (substitutionModel == jc69)
        std::cout << "   * Substitution model                     = " << "JC69" << std::endl;
    else if (substitutionModel == gtr)
        std::cout << "   * Substitution model                     = " << "GTR" << std::endl;
    else
        std::cout << "   * Substitution model                     = " << "Custom" << std::endl;
    if (calculateMarginalLikelihood == true)
        std::cout << "   * Calculate marginal likelihood          = " << "yes" << std::endl;
    else
        std::cout << "   * Calculate marginal likelihood          = " << "no" << std::endl;
    std::cout << "   * Number of MCMC cycles                  = " << numMcmcCycles << std::endl;
    std::cout << "   * Print-to-screen frequency              = " << printFrequency << std::endl;
    std::cout << "   * Chain sample frequency                 = " << sampleFrequency << std::endl;
    std::cout << "   * Tree length prior parameter            = " << inverseTreeLength << std::endl;
    std::cout << std::endl;
}

void UserSettings::usage(void) {

    std::cout << "   Usage" << std::endl;
    std::cout << "   * -d -- File with initial alignments of words" << std::endl;
    std::cout << "   * -m -- Substitution model (jc69/gtr/custom)" << std::endl;
    std::cout << "   * -z -- Calculate marginal likelihood (no/yes)" << std::endl;
    std::cout << "   * -n -- Number of MCMC cycles" << std::endl;
    std::cout << "   * -p -- Print-to-screen frequency" << std::endl;
    std::cout << "   * -s -- Chain sample frequency" << std::endl;
    std::cout << "   * -l -- Inverse of the tree length prior" << std::endl;
    std::cout << std::endl;
}
