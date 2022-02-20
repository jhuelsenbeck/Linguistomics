#ifndef MatrixMath_hpp
#define MatrixMath_hpp

#include "Container.hpp"

#include "Container.hpp"


class MathCache {
    public:
        MathCache(int numstates);
        ~MathCache();

        DoubleMatrix scratch1;
        DoubleMatrix scratch2;
        double* scratchVec;
};

namespace  MatrixMath {

    void    backSubstitutionRow(DoubleMatrix* U, double* b);
    void    computeLandU(DoubleMatrix* A, DoubleMatrix* L, DoubleMatrix* U);
    void    forwardSubstitutionRow(DoubleMatrix* L, double* b);
    void    gaussianElimination(DoubleMatrix* A, DoubleMatrix* B, DoubleMatrix* X, DoubleMatrix* L, DoubleMatrix* U, double* b);
};

#endif
