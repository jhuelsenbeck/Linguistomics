#ifndef ParameterIndelGammaShape_hpp
#define ParameterIndelGammaShape_hpp


#include <string>
#include <vector>
#include "Parameter.hpp"
class Model;
class RandomVariable;



class ParameterIndelGammaShape : public Parameter {

    public:
                                ParameterIndelGammaShape(void) = delete;
                                ParameterIndelGammaShape(const ParameterIndelGammaShape& p) = delete;
                                ParameterIndelGammaShape(RandomVariable* r, Model* m, std::string n, double ep, int nc);
        void                    accept(void);
        std::string             getHeader(void);
        std::string             getJsonString(void);
        std::string             getString(void);
        double                  lnPriorProbability(void);
        std::vector<double>&    getRates(void) { return rates[0]; }
        void                    print(void);
        void                    reject(void);
        double                  update(void);
    
    protected:
        double                  alpha[2];
        double                  expPriorVal;
        int                     numCategories;
        std::vector<double>     rates[2];
};


#endif
