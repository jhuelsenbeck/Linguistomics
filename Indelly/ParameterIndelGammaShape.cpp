#include <cmath>
#include <iomanip>
#include <iostream>
#include "ParameterIndelGammaShape.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "TransitionProbabilities.hpp"



ParameterIndelGammaShape::ParameterIndelGammaShape(RandomVariable* r, Model* m, std::string n, double ep, int nc) : Parameter(r, m, n) {

    std::cout << "   * Setting up gamma shape parameter for insertion/deletion rates " << std::endl;

    updateChangesRateMatrix = false;
    
    expPriorVal = ep;
    numCategories = nc;
    
    alpha[0] = Probability::Exponential::rv(rv, expPriorVal);
    alpha[1] = alpha[0];
    
    rates[0].resize(numCategories);
    rates[1].resize(numCategories);
    Probability::Gamma::discretization(rates[0], alpha[0], alpha[0], numCategories, false);

    rates[1] = rates[0];
}

void ParameterIndelGammaShape::accept(void) {

    alpha[1] = alpha[0];
    rates[1] = rates[0];
}

std::string ParameterIndelGammaShape::getHeader(void) {

    std::string str = "IndelAlpha";
    str += '\t';
    return str;
}

std::string ParameterIndelGammaShape::getString(void) {

    std::string str = std::to_string(alpha[0]);
    str += '\t';
    return str;
}

double ParameterIndelGammaShape::lnPriorProbability(void) {

    return log(expPriorVal) - expPriorVal * alpha[0];
}

void ParameterIndelGammaShape::print(void) {

    std::cout << std::fixed << std::setprecision(6);
    std::cout << alpha[0] << " " << std::endl;
}

void ParameterIndelGammaShape::reject(void) {

    alpha[0] = alpha[1];
    rates[0] = rates[1];
}

double ParameterIndelGammaShape::update(void) {

    lastUpdateType = "indel shape parameter";
    
    double tuning = log(4.0);

    double newAlpha = alpha[0] * exp(tuning*(rv->uniformRv()-0.5));
    double lnP = log(newAlpha) - log(alpha[0]);
    
    alpha[0] = newAlpha;
    Probability::Gamma::discretization(rates[0], alpha[0], alpha[0], numCategories, false);

    updateChangesTransitionProbabilities = true;
    TransitionProbabilities& tip = TransitionProbabilities::transitionProbabilties();
    tip.flipActive();
    tip.setNeedsUpdate(true);
    tip.setTransitionProbabilities();

    return lnP;
}
