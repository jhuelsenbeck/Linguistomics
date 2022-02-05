#include <iomanip>
#include <iostream>
#include "Model.hpp"
#include "Msg.hpp"
#include "ParameterEquilibirumFrequencies.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "RateMatrix.hpp"
#include "TransitionProbabilities.hpp"

double ParameterEquilibirumFrequencies::minVal = 0.001;



ParameterEquilibirumFrequencies::ParameterEquilibirumFrequencies(RandomVariable* r, Model* m, std::string n, int ns) : Parameter(r, m, n) {

    std::cout << "   * Setting up equilibrium frequencies parameter " << std::endl;

    updateChangesRateMatrix = true;

    numStates = ns;
    freqs[0].resize(numStates);
    freqs[1].resize(numStates);
    alpha.resize(numStates);
    for (int i=0; i<numStates; i++)
        alpha[i] = 1.0;
    Probability::Dirichlet::rv(rv, alpha, freqs[0]);
    for (int i=0; i<numStates; i++)
        freqs[0][i] = 1.0 / numStates;
    normalize(freqs[0], minVal);
    freqs[1] = freqs[0];
}

ParameterEquilibirumFrequencies::~ParameterEquilibirumFrequencies(void) {

}

void ParameterEquilibirumFrequencies::accept(void) {

    freqs[1] = freqs[0];
}

std::string ParameterEquilibirumFrequencies::getJsonString(void) {

    std::string str = "";
    for (int i=0; i<numStates; i++)
        {
        str += std::to_string(freqs[0][i]);
        str += '\t';
        }
    return str;
}

std::string ParameterEquilibirumFrequencies::getHeader(void) {

    std::string str = "";
    for (int i=0; i<numStates; i++)
        {
        str += "F[";
        str += std::to_string(i);
        str += "]";
        str += '\t';
        }
    return str;
}

std::string ParameterEquilibirumFrequencies::getString(void) {

    std::string str = "";
    for (int i=0; i<numStates; i++)
        {
        str += std::to_string(freqs[0][i]);
        str += '\t';
        }
    return str;
}

double ParameterEquilibirumFrequencies::lnPriorProbability(void) {

    return Probability::Helper::lnGamma(numStates-1);
}

void ParameterEquilibirumFrequencies::print(void) {

    std::cout << "[ ";
    std::cout << std::fixed << std::setprecision(6);
    for (int i=0; i<numStates; i++)
        {
        std::cout << freqs[0][i] << " ";
        }
    std::cout << "]" << std::endl;
}

void ParameterEquilibirumFrequencies::normalize(std::vector<double>& vec, double min) {

    // find entries with values that are too small
    int numTooSmall = 0;
    double sum = 0.0;
    for (int i=0; i<vec.size(); i++)
        {
        if (vec[i] < min)
            numTooSmall++;
        else
            sum += vec[i];
        }
        
    double factor = (1.0 - numTooSmall * min) / sum;
    for (int i=0; i<vec.size(); i++)
        {
        if (vec[i] < min)
            vec[i] = min;
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

void ParameterEquilibirumFrequencies::normalize(double* vec, double min, int n) {

    // find entries with values that are too small
    int numTooSmall = 0;
    double sum = 0.0;
    for (int i=0; i<n; i++)
        {
        if (vec[i] < min)
            numTooSmall++;
        else
            sum += vec[i];
        }
        
    double factor = (1.0 - numTooSmall * min) / sum;
    for (int i=0; i<n; i++)
        {
        if (vec[i] < min)
            vec[i] = min;
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

std::vector<int> ParameterEquilibirumFrequencies::randomlyChooseIndices(int k, int n) {

    std::vector<int> possibleIndices(n);
    for (int i=0; i<n; i++)
        possibleIndices[i] = i;
        
    std::vector<int> indices(k);
    for (int i=0; i<k; i++)
        {
        int pos = rv->uniformRvInt(n-i);
        indices[i] = possibleIndices[pos];
        possibleIndices[pos] = possibleIndices[n-i-1];
        }
    
    return indices;
}

void ParameterEquilibirumFrequencies::reject(void) {

    freqs[0] = freqs[1];
    modelPtr->flipActiveLikelihood();

}

double ParameterEquilibirumFrequencies::update(void) {

    lastUpdateType = "equilibrium frequencies";

    int k = 1;
    
    double lnP = 0.0;
    if (k == 1)
        {
        double alpha0 = 100.0;
        
        int indexToUpdate = rv->uniformRvInt(numStates);

        oldValues.resize(2, 0.0);
        newValues.resize(2, 0.0);
        alphaForward.resize(2, 0.0);
        alphaReverse.resize(2, 0.0);

        oldValues[0] = freqs[0][indexToUpdate];
        oldValues[1] = 1.0 - oldValues[0];
        alphaForward[0] = oldValues[0] * alpha0;
        alphaForward[1] = oldValues[1] * alpha0;
        
        Probability::Dirichlet::rv(rv, alphaForward, newValues);
        normalize(newValues, minVal);
        
        alphaReverse[0] = newValues[0] * alpha0;
        alphaReverse[1] = newValues[1] * alpha0;
        
        double factor = newValues[1] / oldValues[1];
        
        for (int i=0; i<numStates; i++)
            freqs[0][i] *= factor;
        freqs[0][indexToUpdate] = newValues[0];

        lnP = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
        lnP += (numStates - 2) * log(factor); // Jacobian
        }
    else if (k > 1 && k < numStates)
        {
        double alpha0 = 1000.0;
        std::vector<int> indicesToUpdate = randomlyChooseIndices(k, numStates);
        std::map<size_t,size_t> mapper;
        for (size_t i=0; i<indicesToUpdate.size(); i++)
            mapper.insert( std::make_pair(indicesToUpdate[i], i) );
            
        oldValues.resize(k+1, 0.0);
        newValues.resize(k+1, 0.0);
        alphaForward.resize(k+1, 0.0);
        alphaReverse.resize(k+1, 0.0);
        
        for (size_t i=0; i<numStates; i++)
            {
            std::map<size_t,size_t>::iterator it = mapper.find(i);
            if (it != mapper.end())
                oldValues[it->second] += freqs[0][it->first];
            else
                oldValues[k] += freqs[0][i];
            }
        
        for (size_t i=0; i<k+1; i++)
            alphaForward[i] = oldValues[i] * alpha0;;
        
        // draw a new value for the reduced vector
        Probability::Dirichlet::rv(rv, alphaForward, newValues);
        normalize(newValues, minVal);
        
        // fill in the Dirichlet parameters for the reverse probability calculations
        for (size_t i=0; i<k+1; i++)
            alphaReverse[i] = newValues[i] * alpha0;
        
        // fill in the full vector
        double factor = newValues[k] / oldValues[k];
        for (size_t i=0; i<numStates; i++)
            {
            std::map<size_t,size_t>::iterator it = mapper.find(i);
            if (it != mapper.end())
                freqs[0][i] = newValues[it->second];
            else
                freqs[0][i] *= factor;
            }
        
        // Hastings ratio
        lnP = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
        lnP += (numStates - k - 1) * log(factor); // Jacobian
        }
    else
        {
        // update all of the rates
        double alpha0 = 100000.0;
        oldValues.resize(numStates, 0.0);
        newValues.resize(numStates, 0.0);
        alphaForward.resize(numStates, 0.0);
        alphaReverse.resize(numStates, 0.0);

        std::vector<double> alphaForward(numStates);
        for (int i=0; i<numStates; i++)
            {
            oldValues[i] = freqs[0][i];
            alphaForward[i] = oldValues[i] * alpha0;
            }
        
        std::vector<double> newValues(numStates);
        Probability::Dirichlet::rv(rv, alphaForward, newValues);
        normalize(newValues, minVal);
        
        for (int i=0; i<numStates; i++)
            alphaReverse[i] = newValues[i] * alpha0;

        lnP = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
        
        freqs[0] = newValues;
        }
    
    // update the rate matrix and transition probabilities
    RateMatrix& rmat = RateMatrix::rateMatrix();
    rmat.flipActiveValues();
    rmat.updateRateMatrix(modelPtr->getExchangabilityRates(), freqs[0]);

    updateChangesTransitionProbabilities = true;
    TransitionProbabilities& tip = TransitionProbabilities::transitionProbabilties();
    tip.flipActive();
    tip.setNeedsUpdate(true);
    tip.setTransitionProbabilities();

    modelPtr->setUpdateLikelihood();
    modelPtr->flipActiveLikelihood();

    return lnP;
}

double ParameterEquilibirumFrequencies::updateFromPrior(void) {

    lastUpdateType = "random equilibrium frequencies";

    // draw from the prior distribution, which is a flat Dirichlet distribution
    Probability::Dirichlet::rv(rv, alpha, freqs[0]);

    // update the rate matrix and transition probabilities
    RateMatrix& rmat = RateMatrix::rateMatrix();
    rmat.flipActiveValues();
    rmat.updateRateMatrix(modelPtr->getExchangabilityRates(), freqs[0]);

    updateChangesTransitionProbabilities = true;
    TransitionProbabilities& tip = TransitionProbabilities::transitionProbabilties();
    tip.flipActive();
    tip.setNeedsUpdate(true);
    tip.setTransitionProbabilities();

    return 0.0;
}
