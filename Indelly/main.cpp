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

    // create the thread pool
    ThreadPool threadPool;

    // instantiate the random number generator
    RandomVariable* rv = NULL;
    if (settings.getSeed() == 0)
        rv = new RandomVariable(1);
    else
        rv = new RandomVariable(settings.getSeed());
        
    // print the header
    printHeader();

    // set up the phylogenetic model
    Model** model = new Model*[settings.getNumChains()];
    for (int i=0; i<settings.getNumChains(); i++)
        {
        model[i] = new Model(rv, threadPool);
        model[i]->setIndex(i); 
        }
    
    // run the Markov chain Monte Carlo algorithm
    Mcmc chain(model, rv);
    chain.run();
    
    delete rv;
    for (int i=0; i<settings.getNumChains(); i++)
        delete model[i];
    delete [] model;
    
    return 0;
}

