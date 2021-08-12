//
//  ParameterIndelRates.hpp
//  Indelly
//
//  Created by John Huelsenbeck on 7/30/21.
//  Copyright © 2021 John Huelsenbeck. All rights reserved.
//

#ifndef ParameterIndelRates_hpp
#define ParameterIndelRates_hpp

#include <vector>
#include "Parameter.hpp"
class Model;



class ParameterIndelRates : public Parameter {

    public:
                                        ParameterIndelRates(void) = delete;
                                        ParameterIndelRates(const ParameterIndelRates& pa) = delete;
                                        ParameterIndelRates(RandomVariable* r, Model* m, std::string n, double slen, double insLam, double delLam);
                                       ~ParameterIndelRates(void);
        void                            accept(void);
        double                          getDeletionRate(void);
        double                          getExpectedSequenceLength(void);
        std::string                     getHeader(void);
        double                          getInsertionRate(void);
        std::string                     getString(void);
        double                          lnPriorProbability(void);
        void                            print(void);
        void                            reject(void);
        double                          update(void);
        
    protected:
        double                          expectedEpsilon(double slen);
        double                          updateFromPrior(void);
        double                          insertionLambda;
        double                          deletionLambda;
        double                          expEpsilon;
        std::vector<double>             epsilon[2];
        double                          rhoExpParm;
        double                          rho[2];
    
};

#endif /* ParameterIndelRates_hpp */
