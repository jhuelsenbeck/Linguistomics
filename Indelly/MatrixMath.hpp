#ifndef MatrixMath_hpp
#define MatrixMath_hpp

#include "Container.hpp"

#include "Container.hpp"


class MathCache {
    public:
                     MathCache(int numstates);
                     ~MathCache();
        void         backSubstitutionRow(DoubleMatrix& U, double* b);
        void         computeLandU(DoubleMatrix& A, DoubleMatrix& L, DoubleMatrix& U);
        void         forwardSubstitutionRow(DoubleMatrix& L, double* b);
        void         gaussianElimination(DoubleMatrix& A, DoubleMatrix& B, DoubleMatrix& X, DoubleMatrix& L);
        void         multiply(DoubleMatrix& A, DoubleMatrix& B);


        DoubleMatrix scratch1;
        DoubleMatrix scratch2;
        DoubleMatrix scratchmult;
        double*      scratchVec;

};


#endif
