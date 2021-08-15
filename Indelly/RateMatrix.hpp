#ifndef RateMatrix_hpp
#define RateMatrix_hpp

#include <vector>
#include "EigenSystem.hpp"



class RateMatrix {

    public:
        static RateMatrix& rateMatrix(void) {
                                static RateMatrix m;
                                return m;
                            }
        StateMatrix_t&      getRateMatrix(void) { return Q[activeMatrix]; }
        void                flipActiveValues(void);
        void                initialize(int d, bool useEigens);
        void                updateRateMatrix(std::vector<double>& rates, std::vector<double>& f);
        
    private:
                            RateMatrix(void);
                            RateMatrix(const RateMatrix&) = delete;
        int                 findPadeQValue(const double tol);
        void                setPadeTolerance(const double tol);
        int                 numStates;
        int                 activeMatrix;
        bool                useEigenSystem;
        int                 padeQValue;
        double              padeTolerance;
        StateMatrix_t       Q[2];
        std::vector<double> equilibriumFrequencies[2];
        bool                isInitialized;
};

#endif
