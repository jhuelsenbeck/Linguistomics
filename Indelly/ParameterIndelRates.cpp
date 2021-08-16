#include <cmath>
#include <iomanip>
#include <iostream>
#include "ParameterIndelRates.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "TransitionProbabilities.hpp"



ParameterIndelRates::ParameterIndelRates(RandomVariable* r, Model* m, std::string n, double slen, double insLam, double delLam) : Parameter(r, m, n) {
    
    std::cout << "   * Setting up insertion/deletion rates parameter " << std::endl;

    updateChangesRateMatrix = false;

    insertionLambda = insLam;
    deletionLambda  = delLam;
    
    double expEpsilon = expectedEpsilon(slen); // epsilon = insertionRate / deletionRate
    
    epsilon[0].resize(2);
    epsilon[1].resize(2);
    //std::cout << "Expected Epsilon = " << expEpsilon << std::endl;
    std::vector<double> alpha(2);
    alpha[0] = 100.0 * expEpsilon;
    alpha[1] = 190.0 * (1.0 - expEpsilon);
    Probability::Dirichlet::rv(rv, alpha, epsilon[0]);

    rhoExpParm = 100.0;
    rho[0] = Probability::Exponential::rv(rv, rhoExpParm);
    
    epsilon[1] = epsilon[0];
    rho[1] = rho[0];

    //std::cout << "lambda = " << getInsertionRate() << " mu = " << getDeletionRate() << " E(L) = " << getExpectedSequenceLength() << std::endl;
}

ParameterIndelRates::~ParameterIndelRates(void) {

}

void ParameterIndelRates::accept(void) {

    epsilon[1] = epsilon[0];
    rho[1] = rho[0];
}

double ParameterIndelRates::expectedEpsilon(double slen) {

    double eps = 0.01;
    double increment = 0.001;
    bool increase = true;
    int nSwitches = 0;
    bool stopLoop = false;
    do {
        double expVal = eps / (1.0 - eps);

        if (increase == true && expVal < slen)
            {
            //std::cout << " 1 " << eps << " " << expVal << " " << slen << std::endl;
            increment *= 0.5;
            increase = false;
            nSwitches++;
            eps += increment;
            }
        else if (increase == true && expVal > slen)
            {
            //std::cout << " 2 " << eps << " " << expVal << " " << slen << std::endl;
            eps -= increment;
            }
        else if (increase == false && expVal < slen)
            {
            //std::cout << " 3 " << eps << " " << expVal << " " << slen << std::endl;
            eps += increment;
            }
        else
            {
            //std::cout << " 4 " << eps << " " << expVal << " " << slen << std::endl;
            increment *= 0.5;
            increase = true;
            nSwitches++;
            eps -= increment;
            }
    if (fabs(expVal - slen) < 0.00001)
        stopLoop = true;
    if (nSwitches > 100)
        stopLoop = true;
        
    } while (stopLoop == false);
    return eps;
}

double ParameterIndelRates::getDeletionRate(void) {

    double mu = rho[0] / (1.0 + epsilon[0][0]);
    return mu;
}

double ParameterIndelRates::getExpectedSequenceLength(void) {

    double lambda = getInsertionRate();
    double mu = getDeletionRate();
    double x = lambda / mu;
    return x / (1.0 - x);
}

std::string ParameterIndelRates::getHeader(void) {

    std::string str = "Insertion";
    str += '\t';
    str += "Deletion";
    str += '\t';
    return str;
}

double ParameterIndelRates::getInsertionRate(void) {

    double lambda = rho[0] - getDeletionRate();
    return lambda;
}

std::string ParameterIndelRates::getString(void) {

    std::string str = std::to_string(getInsertionRate());
    str += '\t';
    str += std::to_string(getDeletionRate());
    str += '\t';
    return str;
}

double ParameterIndelRates::lnPriorProbability(void) {

    double lambda = getInsertionRate();
    double mu = getDeletionRate();
    double lnP = log(deletionLambda) + log(insertionLambda + deletionLambda) - deletionLambda * mu - insertionLambda * lambda;
    return lnP;
}

void ParameterIndelRates::print(void) {

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "[ " << getInsertionRate() << " " << getDeletionRate() << " ]" << std::endl;
}

void ParameterIndelRates::reject(void) {

    epsilon[0] = epsilon[1];
    rho[0] = rho[1];
}

double ParameterIndelRates::update(void) {    
    
    double lnProposalProbability = 0.0;
    double u = rv->uniformRv();
    if (u < 0.5)
        {
        // update epsilon
        lastUpdateType = "epsilon";
        double alpha0 = 10.0;
        std::vector<double> forwardAlpha(2);
        forwardAlpha[0] = alpha0 * epsilon[0][0];
        forwardAlpha[1] = alpha0 * epsilon[0][1];
        std::vector<double> newEpsilon(2);
        
        do {
            Probability::Dirichlet::rv(rv, forwardAlpha, newEpsilon);
            } while(newEpsilon[0] > 0.99 || newEpsilon[0] < 0.01);
        
        std::vector<double> reverseAlpha(2);
        reverseAlpha[0] = alpha0 * newEpsilon[0];
        reverseAlpha[1] = alpha0 * newEpsilon[1];
        
        lnProposalProbability = Probability::Dirichlet::lnPdf(reverseAlpha, epsilon[0]) - Probability::Dirichlet::lnPdf(forwardAlpha, newEpsilon);
        epsilon[0] = newEpsilon;
        }
    else
        {
        // update rho
        lastUpdateType = "rho";
        double tuning = log(4.0);
        double newRho = rho[0] * exp(tuning*(rv->uniformRv()-0.5));
        lnProposalProbability = log(newRho) - log(rho[0]);
        rho[0] = newRho;
        }
        
//    std::cout << "1. rho = " << rho[0] << " epsilon = " << epsilon[0][0] << std::endl;
//    std::cout << "2. rho = " << getInsertionRate() + getDeletionRate() << " epsilon = " << getInsertionRate() / getDeletionRate() << std::endl;

    // set flags indicating the transition probabilities are not affected
    TransitionProbabilities& tip = TransitionProbabilities::transitionProbabilties();
    tip.setNeedsUpdate(false);
    updateChangesTransitionProbabilities = false;

    return lnProposalProbability;
}

double ParameterIndelRates::updateFromPrior(void) {

    lastUpdateType = "random indel rates";

    double newLambda = 0.0, newMu = 0.0;
    do
        {
        newLambda = Probability::Exponential::rv(rv, insertionLambda);
        newMu = Probability::Exponential::rv(rv, deletionLambda);
        } while (newLambda > newMu);
        
    double lnP  = log(deletionLambda) + log(insertionLambda + deletionLambda) - deletionLambda * newMu - insertionLambda * newLambda;
           lnP += log(deletionLambda) + log(insertionLambda + deletionLambda) - deletionLambda * getDeletionRate() - insertionLambda * getInsertionRate();
           
    epsilon[0][0] = newLambda / newMu;
    epsilon[0][1] = 1.0 - epsilon[0][0];
    rho[0] = newLambda + newMu;

    return lnP;
}
