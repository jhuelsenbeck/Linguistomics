#include <cmath>
#include <iomanip>
#include <iostream>
#include "ParameterRatesGammaShape.hpp"
#include "RandomVariable.hpp"
#include "TransitionProbabilities.hpp"



ParameterRatesGammaShape::ParameterRatesGammaShape(RandomVariable* r, Model* m, std::string n, double ep, int nc) : Parameter(r, m, n) {

    std::cout << "   * Setting up gamma shape parameter for site rates " << std::endl;

    updateChangesEigens = false;
    
    expPriorVal = ep;
    numCategories = nc;
    
    alpha[0] = rv->exponentialRv(expPriorVal);
    alpha[1] = alpha[0];
    
    rates[0].resize(numCategories);
    rates[1].resize(numCategories);
    rv->discretizeGamma(rates[0], alpha[0], alpha[0], numCategories, false);
    rates[1] = rates[0];
}

void ParameterRatesGammaShape::accept(void) {

    alpha[1] = alpha[0];
    rates[1] = rates[0];
}

std::string ParameterRatesGammaShape::getHeader(void) {

    std::string str = "RatesAlpha";
    str += '\t';
    return str;
}

std::string ParameterRatesGammaShape::getString(void) {

    std::string str = std::to_string(alpha[0]);
    str += '\t';
    return str;
}

double ParameterRatesGammaShape::lnPriorProbability(void) {

    return log(expPriorVal) - expPriorVal * alpha[0];
}

void ParameterRatesGammaShape::print(void) {

    std::cout << std::fixed << std::setprecision(6);
    std::cout << alpha[0] << " " << std::endl;
}

void ParameterRatesGammaShape::reject(void) {

    alpha[0] = alpha[1];
    rates[0] = rates[1];
}

double ParameterRatesGammaShape::update(void) {

    lastUpdateType = "gamma shape parameter";
    
    double tuning = log(4.0);

    double newAlpha = alpha[0] * exp(tuning*(rv->uniformRv()-0.5));
    double lnP = log(newAlpha) - log(alpha[0]);
    
    alpha[0] = newAlpha;
    rv->discretizeGamma(rates[0], alpha[0], alpha[0], numCategories, false);

    updateChangesTransitionProbabilities = true;
    TransitionProbabilities& tip = TransitionProbabilities::transitionProbabilties();
    tip.flipActive();
    tip.setNeedsUpdate(true);
    tip.setTransitionProbabilities();

    return lnP;
}
