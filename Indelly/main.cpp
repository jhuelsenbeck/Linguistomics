#include <iostream>
#include "Mcmc.hpp"
#include "Model.hpp"
#include "RandomVariable.hpp"
#include "UserSettings.hpp"
#include "Threads.hpp"

void printHeader(int count);

#include "Container.hpp"


int main(int argc, char* argv[]) {

    MyDoubleMatrix x(10,10);
    for (int i=0; i<10; i++)
        for (int j=0; j<10; j++)
            x(i, j) = (double)j;
            
    x.print("x");
    
    MyDoubleMatrix a(5,5);
    MyDoubleMatrix b;
    b.print("b");
    b = a;
    
    a.print("a");
    b.print("b");
    
    
    return 0;

    // create the thread pool
    ThreadPool pool;

    // print the header
    printHeader(pool.ThreadCount);

    // read the user settings from the command-line arguments
    UserSettings& settings = UserSettings::userSettings();
    settings.readCommandLineArguments(argc, argv);
    settings.print();
        
    // instantiate the random number generator
    RandomVariable rv;
        
    // set up the phylogenetic model
    Model model(&rv, pool);
    
    // run the Markov chain Monte Carlo algorithm
    Mcmc chain(&model, &rv);
    chain.run();
    
    return 0;
}

void printHeader(int count) {

    std::cout << "   TongueTwister 1.0" << std::endl;
    std::cout << "   * John P. Huelsenbeck (University of California, Berkeley)" << std::endl;
    std::cout << "   * Shawn McCreight (University of California, Berkeley)" << std::endl;
    std::cout << "   * David Goldstein (University of California, Los Angeles)" << std::endl;
    std::cout << "   * Running on " << count << " processors" << std::endl;
    std::cout << std::endl;
}
