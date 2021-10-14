#ifndef RateMatrix_hpp
#define RateMatrix_hpp

#include <vector>
#include "EigenSystem.hpp"



class RateMatrix {

    public:
        static RateMatrix&      rateMatrix(void) {
                                    static RateMatrix m;
                                    return m;
                                }
        std::vector<double>&    getEquilibriumFrequencies(void) { return equilibriumFrequencies[activeMatrix]; }
        StateMatrix_t&          getRateMatrix(void) { return Q[activeMatrix]; }
        const StateMatrix_t&    getRateMatrix(void) const { return Q[activeMatrix]; }
        bool                    getUseEigenSystem(void) { return useEigenSystem; }
        void                    flipActiveValues(void);
        void                    initialize(int d, bool useEigens);
        void                    updateRateMatrix(std::vector<double>& rates, std::vector<double>& f);
        
    private:
                                RateMatrix(void);
                                RateMatrix(const RateMatrix&) = delete;
        int                     numStates;
        int                     activeMatrix;
        bool                    useEigenSystem;
        StateMatrix_t           Q[2];
        std::vector<double>     equilibriumFrequencies[2];
        bool                    isInitialized;
};

#endif
