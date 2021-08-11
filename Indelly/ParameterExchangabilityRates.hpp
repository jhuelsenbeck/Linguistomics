#ifndef ParameterExchangabilityRates_hpp
#define ParameterExchangabilityRates_hpp

#include <vector>
#include "Parameter.hpp"
class Model;



class ParameterExchangabilityRates : public Parameter {

    public:
                                        ParameterExchangabilityRates(void) = delete;
                                        ParameterExchangabilityRates(const ParameterExchangabilityRates& pr) = delete;
                                        ParameterExchangabilityRates(RandomVariable* r, Model* m, std::string n, int ns);
                                        ParameterExchangabilityRates(RandomVariable* r, Model* m, std::string n, int ns, std::vector<std::string> labs);
                                       ~ParameterExchangabilityRates(void);
        void                            accept(void);
        std::string                     getHeader(void);
        std::string                     getString(void);
        std::vector<double>&            getValue(void) { return rates[0]; }
        double                          lnPriorProbability(void);
        void                            print(void);
        void                            reject(void);
        double                          update(void);
        
    private:
        void                            normalize(std::vector<double>& vec, double minVal);
        std::vector<int>                randomlyChooseIndices(int k, int n);
        int                             numStates;
        int                             numRates;
        std::vector<double>             rates[2];
        std::vector<double>             alpha;
        std::vector<std::string>        rateLabels;
        static double                   minVal;
};

#endif
