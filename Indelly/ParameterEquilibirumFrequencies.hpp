#ifndef ParameterEquilibirumFrequencies_hpp
#define ParameterEquilibirumFrequencies_hpp

#include <vector>
#include "Parameter.hpp"
class Model;



class ParameterEquilibirumFrequencies : public Parameter {

    public:
                                        ParameterEquilibirumFrequencies(void) = delete;
                                        ParameterEquilibirumFrequencies(const ParameterEquilibirumFrequencies& pr) = delete;
                                        ParameterEquilibirumFrequencies(RandomVariable* r, Model* m, std::string n, int ns);
                                       ~ParameterEquilibirumFrequencies(void);
        void                            accept(void);
        std::string                     getHeader(void);
        std::vector<double>&            getValue(void) { return freqs[0]; }
        std::string                     getString(void);
        double                          lnPriorProbability(void);
        void                            print(void);
        void                            reject(void);
        double                          update(void);
        
    private:
        void                            normalize(std::vector<double>& vec, double minVal);
        std::vector<int>                randomlyChooseIndices(int k, int n);
        int                             numStates;
        std::vector<double>             freqs[2];
        std::vector<double>             alpha;
};

#endif
