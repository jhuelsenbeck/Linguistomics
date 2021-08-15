#include <iostream>
#include <thread>
#include "Mcmc.hpp"
#include "Model.hpp"
#include "RandomVariable.hpp"
#include "UserSettings.hpp"

void printHeader(void);



int main(int argc, char* argv[]) {

    printHeader();

    // read the user settings from the command-line arguments
    UserSettings& settings = UserSettings::userSettings();
    settings.readCommandLineArguments(argc, argv);
    settings.print();
        
    // instantiate the random number generator
    RandomVariable rv;
        
    // set up the phylogenetic model
    Model model(&rv);
    
    // run the Markov chain Monte Carlo algorithm
    Mcmc chain(&model, &rv);
    chain.run();
    
    return 0;
}

void printHeader(void) {

    std::cout << "   TongueTwister 1.0" << std::endl;
    std::cout << "   * John P. Huelsenbeck (University of California, Berkeley)" << std::endl;
    std::cout << "   * Shawn McCreight (University of California, Berkeley)" << std::endl;
    std::cout << "   * David Goldstein (University of California, Los Angeles)" << std::endl;
    std::cout << "   * Running on " << std::thread::hardware_concurrency() << " processors" << std::endl;
    std::cout << std::endl;
}
