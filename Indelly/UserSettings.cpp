#include <fstream>
#include <iostream>
#include <vector>
#include "json.hpp"
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
    numRateCategories           = 1;
    numIndelCategories          = 1;
    useEigenSystem              = true;
    useOnlyCompleteWords        = false;
}

void UserSettings::readCommandLineArguments(int argc, char* argv[]) {

    // move the command-line arguments to a good-ol' C++ STL vector of strings
    std::vector<std::string> commands;
    executablePath = argv[0];
#   if defined(DEBUG_MODE)
    commands.push_back(executablePath);
    commands.push_back("-d");
    commands.push_back("/Users/johnh/GitHub/Linguistomics/Run/Run1/config.json");
    commands.push_back("-c");
    commands.push_back("no");
    commands.push_back("-o");
    commands.push_back("/Users/johnh/Desktop/Indelly/out/test_custom");
    commands.push_back("-m");
    commands.push_back("custom");
    commands.push_back("-n");
    commands.push_back("100");
    commands.push_back("-p");
    commands.push_back("1");
    commands.push_back("-s");
    commands.push_back("100");
    commands.push_back("-l");
    commands.push_back("0.15");
    commands.push_back("-z");
    commands.push_back("no");
    commands.push_back("-e");
    commands.push_back("no");
#   else
    for (int i=0; i<argc; i++)
        commands.push_back(argv[i]);
#   endif
    
    // check that we have an odd number of commands
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
            else if (arg == "-e")
                {
                if (cmd == "yes")
                    useEigenSystem = true;
                else if (cmd == "no")
                    useEigenSystem = false;
                else
                    Msg::error("Unknon option for calculating matrix exponential");
                }
            else if (arg == "-c")
                {
                if (cmd == "yes")
                    useOnlyCompleteWords = true;
                else if (cmd == "no")
                    useOnlyCompleteWords = false;
                else
                    Msg::error("Unknon option for use of complete words");
                }
            else if (arg == "-g")
                numRateCategories = atoi(cmd.c_str());
            else if (arg == "-i")
                numIndelCategories = atoi(cmd.c_str());
            else
                {
                usage();
                Msg::error("Unknown command \"" + arg + "\"");
                }
            
            arg = "";
            }
        }
        
    // check if the data file contains the commands. If so, we will use those, overwriting what we
    // just parsed (excepting the path/name to the configuration file
    std::string fn = getDataFile();
    std::ifstream ifs(fn);
    nlohmann::json j;
    try
        {
        j = nlohmann::json::parse(ifs);
        }
    catch (nlohmann::json::parse_error& ex)
        {
        Msg::error("Error parsing JSON file at byte " + std::to_string(ex.byte));
        }
    auto it = j.find("McmcSettings");
    if (it != j.end())
        {
        // reading settings from program
        nlohmann::json jsonSettings = j["McmcSettings"];
        // output file (FileOutput)
        // only analyze complete words yes/no (OnlyCompleteWords)
        // substitution model JC69/Custom/GTR (Model)
        // calculate marginal likelihood yes/no (CalcMarginal)
        // use eigen system for matrix exponential yes/no (UseEigenSystem)
        // number of mcmc cycles (NumCycles)
        // print to screen frequency (PrintFreq)
        // sample frequency (SampleFreq)
        // tree-length prior parameter value (TreeLengthPriorVal)
        
        }
}

void UserSettings::print(void) {

    std::cout << "   Settings" << std::endl;
    std::cout << "   * Executable path                         = " << executablePath << std::endl;
    std::cout << "   * File with initial word alignments       = \"" << dataFile << "\"" << std::endl;
    if (useOnlyCompleteWords == true)
        std::cout << "   * Only analyzing completely sampled words = yes" << std::endl;
    else
        std::cout << "   * Only analyzing completely sampled words = no" << std::endl;
    if (substitutionModel == jc69)
        std::cout << "   * Substitution model                      = JC69" << std::endl;
    else if (substitutionModel == gtr)
        std::cout << "   * Substitution model                      = GTR" << std::endl;
    else
        std::cout << "   * Substitution model                      = Custom" << std::endl;
    std::cout << "   * Number of gamma rate categories         = " << numRateCategories << std::endl;
    std::cout << "   * Number of gamma indel categories        = " << numIndelCategories << std::endl;
    if (calculateMarginalLikelihood == true)
        std::cout << "   * Calculate marginal likelihood           = yes" << std::endl;
    else
        std::cout << "   * Calculate marginal likelihood           = no" << std::endl;
    if (useEigenSystem == true)
        std::cout << "   * Use eigen system for matrix exponential = yes" << std::endl;
    else
        std::cout << "   * Use eigen system for matrix exponential = no" << std::endl;
    std::cout << "   * Number of MCMC cycles                   = " << numMcmcCycles << std::endl;
    std::cout << "   * Print-to-screen frequency               = " << printFrequency << std::endl;
    std::cout << "   * Chain sample frequency                  = " << sampleFrequency << std::endl;
    std::cout << "   * Tree length prior parameter             = " << inverseTreeLength << std::endl;
    std::cout << std::endl;
}

void UserSettings::usage(void) {

    std::cout << "   Usage" << std::endl;
    std::cout << "   * -d -- File with initial alignments of words" << std::endl;
    std::cout << "   * -m -- Substitution model (jc69/gtr/custom)" << std::endl;
    std::cout << "   * -c -- Use only completely sampled words (no/yes)" << std::endl;
    std::cout << "   * -g -- Number of gamma rate categories (=1 is no rate variation)" << std::endl;
    std::cout << "   * -i -- Number of gamma indel categories (=1 is no indel rate variation)" << std::endl;
    std::cout << "   * -z -- Calculate marginal likelihood (no/yes)" << std::endl;
    std::cout << "   * -n -- Number of MCMC cycles" << std::endl;
    std::cout << "   * -p -- Print-to-screen frequency" << std::endl;
    std::cout << "   * -s -- Chain sample frequency" << std::endl;
    std::cout << "   * -e -- Use the Eigen system when calculating the matrix exponential (no/yes)" << std::endl;
    std::cout << "   * -l -- Inverse of the tree length prior" << std::endl;
    std::cout << std::endl;
}
