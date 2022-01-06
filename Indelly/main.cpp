#include <iostream>
#include "Mcmc.hpp"
#include "Model.hpp"
#include "RandomVariable.hpp"
#include "UserSettings.hpp"
#include "thread_pool.hpp"

void printHeader(thread_pool* p);



int main(int argc, char* argv[]) {

    // create the thread pool
    thread_pool pool;

    // print hte header
    printHeader(&pool);

    // read the user settings from the command-line arguments
    UserSettings& settings = UserSettings::userSettings();
    settings.readCommandLineArguments(argc, argv);
    settings.print();
        
    // instantiate the random number generator
    RandomVariable rv(1);
        
    // set up the phylogenetic model
    Model model(&rv, &pool);
    
    // run the Markov chain Monte Carlo algorithm
    Mcmc chain(&model, &rv);
    chain.run();
    
    return 0;
}

void printHeader(thread_pool* p) {

    std::cout << "   TongueTwister 1.0" << std::endl;
    std::cout << "   * John P. Huelsenbeck (University of California, Berkeley)" << std::endl;
    std::cout << "   * Shawn McCreight (University of California, Berkeley)" << std::endl;
    std::cout << "   * David Goldstein (University of California, Los Angeles)" << std::endl;
    std::cout << "   * Running on " << p->get_thread_count() << " processors" << std::endl;
    std::cout << std::endl;
}
