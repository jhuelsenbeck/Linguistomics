#ifndef RateMatrix_hpp
#define RateMatrix_hpp

#include <vector>
#include "Container.hpp"



class RateMatrix {

    public:
        static RateMatrix&      rateMatrix(void) {
                                    static RateMatrix m;
                                    return m;
                                }
        std::vector<double>&    getEquilibriumFrequencies(void) { return equilibriumFrequencies[activeMatrix]; }
        DoubleMatrix&           getRateMatrix(void) { return *Q[activeMatrix]; }
        void                    flipActiveValues(void);
        void                    initialize(int d);
        void                    updateRateMatrix(std::vector<double>& rates, std::vector<double>& f);
        
    private:
                                RateMatrix(void);
                                RateMatrix(const RateMatrix&) = delete;
                               ~RateMatrix(void);
        int                     numStates;
        int                     activeMatrix;
        DoubleMatrix*           Q[2];
        std::vector<double>     equilibriumFrequencies[2];
        bool                    isInitialized;
};

#endif
