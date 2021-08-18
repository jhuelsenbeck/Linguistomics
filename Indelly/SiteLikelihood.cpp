#include <iomanip>
#include "SiteLikelihood.hpp"



SiteLikelihood::SiteLikelihood(int nn, int ns) {

    // set dimensions
    numNodes = nn;
    numStates = ns;
    
    // allocate arrays
    zeroH = new double[numStates];
    for (int i=0; i<numStates; i++)
        zeroH[i] = 0.0;
        
    zeroI = new double[numStates + 1];
    for (int i=0; i<numStates+1; i++)
        zeroI[i] = 0.0;

    probsH = new double*[numNodes];
    probsH[0] = new double[numNodes * numStates];
    for (int i=0; i<numNodes; i++)
        probsH[i] = probsH[i-1] + numStates;
    
    probsI = new double*[numNodes];
    probsI[0] = new double[numNodes * (numStates + 1)];
    for (int i=0; i<numNodes; i++)
        probsI[i] = probsI[i-1] + (numStates + 1);
                    
    // set arrays
    for (int i=0; i<numNodes; i++)
        {
        zeroOutH(i);
        zeroOutI(i);
        }
}

SiteLikelihood::~SiteLikelihood(void) {

    delete [] probsH[0];
    delete [] probsH;
    delete [] probsI[0];
    delete [] probsI;
    delete [] zeroH;
    delete [] zeroI;
}

void SiteLikelihood::print(void) {

    std::cout << std::fixed << std::setprecision(1);
    
    std::cout << "probsH vector:" << std::endl;
    for (int i=0; i<numNodes; i++)
        {
        std::cout << std::setw(3) << i << " -- ";
        for (int j=0; j<numStates; j++)
            std::cout << probsH[i][j] << " ";
        std::cout << std::endl;
        }

    std::cout << "probsI vector:" << std::endl;
    for (int i=0; i<numNodes; i++)
        {
        std::cout << std::setw(3) << i << " -- ";
        for (int j=0; j<numStates+1; j++)
            std::cout << probsI[i][j] << " ";
        std::cout << std::endl;
        }
}
