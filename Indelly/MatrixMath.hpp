#ifndef MatrixMath_hpp
#define MatrixMath_hpp

#include "Container.hpp"

namespace  MatrixMath {

    void    backSubstitutionRow(DoubleMatrix* U, double* b);
    void    computeLandU(DoubleMatrix* A, DoubleMatrix* L, DoubleMatrix* U);
    void    divideMatrixByPowerOfTwo(DoubleMatrix* M, int power);
    void    forwardSubstitutionRow(DoubleMatrix* L, double* b);
    void    gaussianElimination(DoubleMatrix* A, DoubleMatrix* B, DoubleMatrix* X, DoubleMatrix* L, DoubleMatrix* U, double* b);
    void    multiplicationByScalar(DoubleMatrix* M, double c);
    void    multiplicationByScalar(DoubleMatrix* M, double c, DoubleMatrix* Res);
    void    multiplyTwoMatrices(const DoubleMatrix* A, const DoubleMatrix* B, DoubleMatrix* C, DoubleMatrix* temp);
};

#endif
