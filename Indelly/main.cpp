#include <iostream>
#include "Mcmc.hpp"
#include "Model.hpp"
#include "RandomVariable.hpp"
#include "UserSettings.hpp"
#include "Threads.hpp"

void printHeader() {

    std::cout << "   TongueTwister 1.0" << std::endl;
    std::cout << "   * John P. Huelsenbeck (University of California, Berkeley)" << std::endl;
    std::cout << "   * Shawn McCreight (University of California, Berkeley)" << std::endl;
    std::cout << "   * David Goldstein (University of California, Los Angeles)" << std::endl;
    std::cout << std::endl;
}


int main(int argc, char* argv[]) {

    // read the user settings from the command-line arguments
    UserSettings& settings = UserSettings::userSettings();
    settings.readCommandLineArguments(argc, argv);
    settings.print();
        
    // instantiate the random number generator
    RandomVariable rv;
        
    // print the header
    printHeader();

    // set up the phylogenetic model
    Model model(&rv);
    
    // run the Markov chain Monte Carlo algorithm
    Mcmc chain(&model, &rv);
    chain.run();
    
    return 0;
}

