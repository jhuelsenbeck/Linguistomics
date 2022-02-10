#ifndef MatrixMath_hpp
#define MatrixMath_hpp

#include "JphMatrix.hpp"

namespace  MatrixMath {

    void    addTwoMatrices(DoubleMatrix* A, DoubleMatrix* B, DoubleMatrix* Res);
    void    backSubstitutionRow(DoubleMatrix* U, double* b);
    void    computeLandU(DoubleMatrix* A, DoubleMatrix* L, DoubleMatrix* U);
    void    divideMatrixByPowerOfTwo(DoubleMatrix* M, int power);
    void    forwardSubstitutionRow(DoubleMatrix* L, double* b);
    void    gaussianElimination(DoubleMatrix* A, DoubleMatrix* B, DoubleMatrix* X, DoubleMatrix* L, DoubleMatrix* U, double* b);
    void    multiplicationByScalar(DoubleMatrix* M, double c);
    void    multiplicationByScalar(DoubleMatrix* M, double c, DoubleMatrix* Res);
    void    setIdentity(DoubleMatrix* M);
    void    multiplyTwoMatrices(DoubleMatrix* A, DoubleMatrix* B, DoubleMatrix* C, DoubleMatrix* temp);
};

#endif
