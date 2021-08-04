#ifndef EigenSystem_H
#define EigenSystem_H

#include <complex>
#include <vector>
#include <Eigen/Core>

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> StateVector_t;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> StateMatrix_t;



class EigenSystem {

    public:
        static EigenSystem&                                     eigenSystem(void)
                                                                    {
                                                                    static EigenSystem singleEigenSystem;
                                                                    return singleEigenSystem;
                                                                    }
        void                                                    calculateEigenSystem(StateMatrix_t& Q);
        std::vector<double>                                     calulateStationaryFrequencies(StateMatrix_t& Q);
        std::vector<double>&                                    getStationaryFrequencies(void) { return pi[activeVals]; }
        void                                                    flipActiveValues(void);
        std::complex<double>*                                   getCijk(void);
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>& getEigenValues(void);
        void                                                    initialize(int ns);
        void                                                    setStationaryFrequencies(std::vector<double> f) { pi[activeVals] = f; }
        void                                                    testStationaryFrequencies(StateMatrix_t& Q);
        void                                                    updateRateMatrix(std::vector<double>& rates, std::vector<double>& f);

    private:
                                                                EigenSystem(void);
                                                                EigenSystem(const EigenSystem&);
                                                                EigenSystem& operator=(const EigenSystem&);
                                                               ~EigenSystem(void);
        double                                                  sumVector(std::vector<double>& v);
        bool                                                    isInitialized;
        int                                                     numStates;
        int                                                     activeVals;
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>  eigenValues[2];
        std::complex<double>*                                   ccIjk[2];
        std::vector<double>                                     pi[2];
};

#endif
