#include <cmath>
#include <iomanip>
#include <iostream>
#include "RateMatrix.hpp"
#include "RateMatrixHelper.hpp"



RateMatrix::RateMatrix(void) {

    numStates = 0;
    useEigenSystem = true;
    activeMatrix = 0;
    isInitialized = false;
}

void RateMatrix::flipActiveValues(void) {

    if (activeMatrix == 0)
        activeMatrix = 1;
    else
        activeMatrix = 0;
}

void RateMatrix::initialize(int d, bool useEigens) {

    if (isInitialized == true)
        return;
        
    numStates = d;
    useEigenSystem = useEigens;
    activeMatrix = 0;
    Q[0].resize(numStates,numStates);
    Q[1].resize(numStates,numStates);
    equilibriumFrequencies[0].resize(numStates);
    equilibriumFrequencies[1].resize(numStates);

    isInitialized = true;
}

void  RateMatrix::updateRateMatrix(std::vector<double>& rates, std::vector<double>& f) {

    // set the stationary frequencies
    equilibriumFrequencies[activeMatrix] = f;
    
    // initialize the rate matrix
    StateMatrix_t& Q = this->Q[activeMatrix];
    
    // fill in off diagonal components of the rate matrix in
    // a model-dependent manner
    if (rates.size() == numStates * (numStates-1) / 2)
        {
        // gtr model
        for (int i=0, k=0; i<numStates; i++)
            {
            for (int j=i+1; j<numStates; j++)
                {
                Q(i,j) = rates[k] * f[j];
                Q(j,i) = rates[k] * f[i];
                k++;
                }
            }
        }
    else
        {
        // custom model
        RateMatrixHelper& helper = RateMatrixHelper::rateMatrixHelper();
        int** map = helper.getMap();
        for (int i=0; i<numStates; i++)
            {
            for (int j=i+1; j<numStates; j++)
                {
                int changeType = map[i][j];
                Q(i,j) = rates[changeType] * f[j];
                Q(j,i) = rates[changeType] * f[i];
                }
            }
        }
                
    // fill in the diagonal elements of the rate matrix
    for (int i=0; i<numStates; i++)
        {
        double sum = 0.0;
        for (int j=0; j<numStates; j++)
            {
            if (i != j)
                sum += Q(i,j);
            }
        Q(i,i) = -sum;
        }
    
    // rescale the rate matrix
    double averageRate = 0.0;
    for (int i=0; i<numStates; i++)
        averageRate += -f[i] * Q(i,i);
    double scaleFactor = 1.0 / averageRate;
    Q *= scaleFactor;
    
    // update the eigen system
    if (useEigenSystem == true)
        {
        EigenSystem& eigs = EigenSystem::eigenSystem();
        eigs.setActiveEigens(activeMatrix);
        eigs.calculateEigenSystem(Q);
        }
        
#   if 0
    std::cout << std::fixed << std::setprecision(5);
    std::cout << Q << std::endl;
    for (int i=0; i<numStates; i++)
        std::cout << equilibriumFrequencies[activeMatrix][i] << " ";
    std::cout << std::endl;
#   endif
}
