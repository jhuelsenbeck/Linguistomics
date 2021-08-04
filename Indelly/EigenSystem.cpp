#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include <iomanip>
#include <iostream>
#include <vector>
#include "EigenSystem.hpp"
#include "Msg.hpp"



EigenSystem::EigenSystem(void) {

    isInitialized = false;
}

EigenSystem::~EigenSystem(void) {

    delete [] ccIjk[0];
    delete [] ccIjk[1];
}

void EigenSystem::calculateEigenSystem(StateMatrix_t& Q) {

    // calculate the Eigenvalues and Eigenvectors and do some precomputation
    Eigen::EigenSolver< StateMatrix_t > eigSolver;
    eigSolver.compute( Q, true );
    eigenValues[activeVals] = eigSolver.eigenvalues();
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> eigenVectors = eigSolver.eigenvectors();
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> inverseEigenVectors = eigenVectors.inverse();

    // calculate cc_ijk
    std::complex<double>* pc = ccIjk[activeVals];
    for (int i=0; i<numStates; i++)
        for (int j=0; j<numStates; j++)
            for (int k=0; k<numStates; k++)
                 *(pc++) = eigenVectors(i,k) * inverseEigenVectors(k,j);
}

std::vector<double> EigenSystem::calulateStationaryFrequencies(StateMatrix_t& Q) {
    
    Eigen::VectorXd f = Q.transpose().fullPivLu().kernel();
    f = f / f.sum();
    
    std::vector<double> stationaryFrequencies(numStates);
    for (int i=0; i<numStates; i++)
        stationaryFrequencies[i] = f(i);
    
    double sum = sumVector(stationaryFrequencies);
    if ( fabs(1.0 - sum) > 0.001 )
        {
        std::cout << "sum = " << sum << std::endl;
        Msg::error("Stationary frequencies don't sum to one");
        }
    for (int i=0; i<numStates; i++)
        {
        if (f[i] < 0.0)
            {
            std::cout << f << std::endl;
            Msg::error("Negative stationary frequency");
            }
        }
    
    return stationaryFrequencies;
}

void EigenSystem::flipActiveValues(void) {

    if (activeVals == 0)
        activeVals = 1;
    else
        activeVals = 0;
}

std::complex<double>* EigenSystem::getCijk(void) {

    return ccIjk[activeVals];
}

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>& EigenSystem::getEigenValues(void) {

    return eigenValues[activeVals];
}

void EigenSystem::initialize(int ns) {

    if (isInitialized == true)
        {
        Msg::warning("EigenSystem is already initialized");
        return;
        }
        
    numStates = ns;
    activeVals = 0;
    ccIjk[0] = new std::complex<double>[numStates*numStates*numStates];
    ccIjk[1] = new std::complex<double>[numStates*numStates*numStates];
    pi[0].resize(numStates);
    pi[1].resize(numStates);

    isInitialized = true;
}

double EigenSystem::sumVector(std::vector<double>& v) {

    double sum = 0.0;
    for (int i=0; i<numStates; i++)
        sum += v[i];
    return sum;
}

void EigenSystem::testStationaryFrequencies(StateMatrix_t& Q) {

}

void  EigenSystem::updateRateMatrix(std::vector<double>& rates, std::vector<double>& f) {

    // initialize the rate matrix
    StateMatrix_t Q(numStates,numStates);
        
    // fill in off diagonal components of rate matrix
    for (int i=0, k=0; i<numStates; i++)
        {
        for (int j=i+1; j<numStates; j++)
            {
            Q(i,j) = rates[k] * f[j];
            Q(j,i) = rates[k] * f[i];
            k++;
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
    
    // calculate the Eigenvalues and Eigenvectors and do some precomputation
    calculateEigenSystem(Q);
    setStationaryFrequencies(f);
    
#   if 0
    std::cout << std::fixed << std::setprecision(5);
    std::cout << Q << std::endl;
    for (int i=0; i<numStates; i++)
        std::cout << sf[i] << " ";
    std::cout << std::endl;
#   endif
}
