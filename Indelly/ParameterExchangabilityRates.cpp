#include <iomanip>
#include <iostream>
#include <map>
#include "EigenSystem.hpp"
#include "Model.hpp"
#include "ParameterExchangabilityRates.hpp"
#include "RandomVariable.hpp"
#include "TransitionProbabilities.hpp"



ParameterExchangabilityRates::ParameterExchangabilityRates(RandomVariable* r, Model* m, std::string n, int ns, std::string s) : Parameter(r, m, n) {

    updateChangesEigens = true;

    numStates = ns;
    states = s;
    numRates = numStates * (numStates-1) / 2;
    rates[0].resize(numRates);
    rates[1].resize(numRates);
    alpha.resize(numRates);
    for (int i=0; i<numRates; i++)
        alpha[i] = 1.0;
    rv->dirichletRv(alpha, rates[0]);
    rates[1] = rates[0];
    
}

ParameterExchangabilityRates::~ParameterExchangabilityRates(void) {

}

void ParameterExchangabilityRates::accept(void) {

    rates[1] = rates[0];
}

std::string ParameterExchangabilityRates::getHeader(void) {

    std::string str = "";
    for (int i=0; i<numStates; i++)
        {
        for (int j=i+1; j<numStates; j++)
            {
            str += "R[";
            str += states[i];
            str += "-";
            str += states[j];
            str += "]";
            str += '\t';
            }
        }
    return str;
}

std::string ParameterExchangabilityRates::getString(void) {

    std::string str = "";
    for (int i=0; i<numRates; i++)
        {
        str += std::to_string(rates[0][i]);
        str += '\t';
        }
    return str;
}

double ParameterExchangabilityRates::lnPriorProbability(void) {

    return rv->lnGamma(numRates-1);
}

void ParameterExchangabilityRates::normalize(std::vector<double>& vec, double minVal) {

    // find entries with values that are too small
    int numTooSmall = 0;
    double sum = 0.0;
    for (int i=0; i<vec.size(); i++)
        {
        if (vec[i] < minVal)
            numTooSmall++;
        else
            sum += vec[i];
        }
        
    double factor = (1.0 - numTooSmall * minVal) / sum;
    for (int i=0; i<vec.size(); i++)
        {
        if (vec[i] < minVal)
            vec[i] = minVal;
        else
            vec[i] *= factor;
        }
        
#   if 0
    sum = 0.0;
    for (int i=0; i<vec.size(); i++)
        sum += vec[i];
    if ( fabs(1.0 - sum) > 0.000001)
        std::cout << "Problem normalizing vector " << std::fixed << std::setprecision(20) << sum << std::endl;
#   endif
}

std::vector<int> ParameterExchangabilityRates::randomlyChooseIndices(int k, int n) {

    std::vector<int> possibleIndices(n);
    for (int i=0; i<n; i++)
        possibleIndices[i] = i;
        
    std::vector<int> indices(k);
    for (int i=0; i<k; i++)
        {
        int pos = (int)(rv->uniformRv()*(n-i));
        indices[i] = possibleIndices[pos];
        possibleIndices[pos] = possibleIndices[n-i-1];
        }
    
    return indices;
}

void ParameterExchangabilityRates::reject(void) {

    rates[0] = rates[1];
}

double ParameterExchangabilityRates::update(void) {
    
    int k = 1;
    
    double lnP = 0.0;
    if (k == 1)
        {
        double alpha0 = 1000.0;
        
        int indexToUpdate = (int)(rv->uniformRv()*numRates);

        std::vector<double> oldValues(2, 0.0);
        std::vector<double> newValues(2, 0.0);
        std::vector<double> alphaForward(2, 0.0);
        std::vector<double> alphaReverse(2, 0.0);

        oldValues[0] = rates[0][indexToUpdate];
        oldValues[1] = 1.0 - oldValues[0];
        alphaForward[0] = oldValues[0] * alpha0;
        alphaForward[1] = oldValues[1] * alpha0;
        
        rv->dirichletRv(alphaForward, newValues);
        normalize(newValues, 0.00000001);
        
        alphaReverse[0] = newValues[0] * alpha0;
        alphaReverse[1] = newValues[1] * alpha0;
        
        double factor = newValues[1] / oldValues[1];
        
        for (int i=0; i<numRates; i++)
            rates[0][i] *= factor;
        rates[0][indexToUpdate] = newValues[0];

        lnP = rv->lnDirichletPdf(alphaReverse, oldValues) - rv->lnDirichletPdf(alphaForward, newValues);
        lnP += (numRates - 2) * log(factor); // Jacobian
        }
    else if (k > 1 && k < numRates)
        {
        double alpha0 = 1000.0;
        std::vector<int> indicesToUpdate = randomlyChooseIndices(k, numRates);
        std::map<size_t,size_t> mapper;
        for (size_t i=0; i<indicesToUpdate.size(); i++)
            mapper.insert( std::make_pair(indicesToUpdate[i], i) );
            
        std::vector<double> oldValues(k+1, 0.0);
        std::vector<double> newValues(k+1, 0.0);
        std::vector<double> alphaForward(k+1, 0.0);
        std::vector<double> alphaReverse(k+1, 0.0);
        
        for (size_t i=0; i<numRates; i++)
            {
            std::map<size_t,size_t>::iterator it = mapper.find(i);
            if (it != mapper.end())
                oldValues[it->second] += rates[0][it->first];
            else
                oldValues[k] += rates[0][i];
            }
        
        for (size_t i=0; i<k+1; i++)
            alphaForward[i] = oldValues[i] * alpha0;;
        
        // draw a new value for the reduced vector
        rv->dirichletRv(alphaForward, newValues);
        normalize(newValues, 0.00000001);
        
        // fill in the Dirichlet parameters for the reverse probability calculations
        for (size_t i=0; i<k+1; i++)
            alphaReverse[i] = newValues[i] * alpha0;
        
        // fill in the full vector
        double factor = newValues[k] / oldValues[k];
        for (size_t i=0; i<numRates; i++)
            {
            std::map<size_t,size_t>::iterator it = mapper.find(i);
            if (it != mapper.end())
                rates[0][i] = newValues[it->second];
            else
                rates[0][i] *= factor;
            }
        
        // Hastings ratio
        lnP = rv->lnDirichletPdf(alphaReverse, oldValues) - rv->lnDirichletPdf(alphaForward, newValues);
        lnP += (numRates - k - 1) * log(factor); // Jacobian
//        std::cout << "lnP 2 = " << lnP << std::endl;
//        std::cout << std::fixed << std::setprecision(10);
//        std::cout << "factor = " << factor << std::endl;
//        std::cout << "log(factor) = " << log(factor) << std::endl;
//        for (int i=0; i<k+1; i++)
//            std::cout << alphaForward[i] << " " << oldValues[i] << " -> " << newValues[i] << " " << alphaReverse[i] << std::endl;
        }
    else
        {
        double alpha0 = 100000.0;
        // update all of the rates
        std::vector<double>& oldValues = rates[0];
        std::vector<double> alphaForward(numRates);
        for (int i=0; i<numRates; i++)
            alphaForward[i] = oldValues[i] * alpha0;
        
        std::vector<double> newValues(numRates);
        rv->dirichletRv(alphaForward, newValues);
        normalize(newValues, 0.000001);
        
        std::vector<double> alphaReverse(numRates);
        for (int i=0; i<numRates; i++)
            alphaReverse[i] = newValues[i] * alpha0;

        lnP = rv->lnDirichletPdf(alphaReverse, oldValues) - rv->lnDirichletPdf(alphaForward, newValues);
        
        rates[0] = newValues;
        }
    
    // update the eigen system and transition probabilities
    EigenSystem& eigs = EigenSystem::eigenSystem();
    eigs.flipActiveValues();
    eigs.updateRateMatrix(rates[0], modelPtr->getEquilibriumFrequencies());

    updateChangesTransitionProbabilities = true;
    TransitionProbabilities& tip = TransitionProbabilities::transitionProbabilties();
    tip.flipActive();
    tip.setNeedsUpdate(true);
    tip.setTransitionProbabilities();
    
    // and it also changes the eigen system

    return lnP;
}
