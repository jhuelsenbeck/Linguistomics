#ifndef RandomVariable_hpp
#define RandomVariable_hpp

#include <cstdint>
#include <vector>

enum {

    N = 624,
    M = 397,
    R = 31,
    A = 0x9908B0DF,
    F = 1812433253,
    U = 11,
    S = 7,
    B = 0x9D2C5680,
    T = 15,
    C = 0xEFC60000,
    L = 18,
    MASK_LOWER = (1ull << R) - 1,
    MASK_UPPER = (1ull << R)
};

class RandomVariable {

    public:
                    RandomVariable(void);
                    RandomVariable(uint32_t seed);
        double      betaRv(double a, double b);
        double      betaPdf(double a, double b, double x);
        double      lnBetaPdf(double a, double b, double x);
        double      betaCdf(double a, double b, double x);
        double      betaQuantile(double a, double b, double p);
        void        dirichletRv(const std::vector<double> &a, std::vector<double> &z);
        void        discretizeGamma(std::vector<double>& catRate, double a, double b, int nCats, bool median);
        double      exponentialRv(double lambda);
        double      gammaRv(double a, double b);
        bool        isNan(double x);
        double      lnDirichletPdf(const std::vector<double> &a, const std::vector<double> &z);
        double      lnExponentialPdf(double lambda, double x);
        double      lnGamma(double a);
        double      lnGammaPdf(double a, double b, double x);
        double      uniformRv(void);

    protected:
        double      beta(double a, double b);
        void        buildRandomTree(void);
        double      chebyshev_eval(double x, const double *a, const int n);
        double      chiSquareQuantile(double prob, double v);
        uint32_t    extractU32(void);
        double      gamma(double x);
        double      incompleteBeta(double a, double b, double x);
        double      incompleteGamma (double x, double alpha, double LnGamma_alpha);
        void        initialize(uint32_t seed);
        bool        isFinite(double x);
        double      lnBeta(double a, double b);
        double      lnGammacor(double x);
        double      log1p(double x);
        double      pointNormal(double prob);
        double      rndGamma(double s);
        double      rndGamma1(double s);
        double      rndGamma2(double s);
        double      stirlerr(double n);
        void        twist(void);
        uint16_t    index;
        uint32_t    mt[N];
};

#endif
