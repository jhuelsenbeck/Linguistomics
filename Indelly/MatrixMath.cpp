#include "MatrixMath.hpp"
#include "Msg.hpp"

#define ROW_MAJOR_ORDER


void MatrixMath::addTwoMatrices(DoubleMatrix* A, DoubleMatrix* B, DoubleMatrix* Res) {

#   if defined(ROW_MAJOR_ORDER)

    if (A->size() != B->size() || A->size() != Res->size())
        Msg::error("Problems with matrix dimensions in addTwoMatrices");
        
    for (auto aPtr=A->begin(), bPtr=B->begin(), resPtr=Res->begin(); aPtr != A->end(); aPtr++, bPtr++, resPtr++)
        *resPtr = *aPtr + *bPtr;
        
#   else

    int nr = (int)A->getNumRows();
    int nc = (int)A->getNumCols();
    
    if (nr != B->getNumRows() || nc != B->getNumCols())
        Msg::error("Problems with matrix dimensions in addTwoMatrices");
        
    for (int i=0; i<nr; i++)
        for (int j=0; j<nc; j++)
            (*Res)[i][j] = (*A)[i][j] + (*B)[i][j];
            
#   endif
}

void MatrixMath::backSubstitutionRow(DoubleMatrix* U, double* b) {

    int n = (int)U->getNumRows();
    
    b[n-1] /= (*U)(n-1,n-1);
    for (int i=n-2; i>=0; i--)
        {
        double dotProduct = 0.0;
        for (int j=i+1; j<n; j++)
            dotProduct += (*U)(i,j) * b[j];
        b[i] = (b[i] - dotProduct) / (*U)(i,i);
        }
}

void MatrixMath::computeLandU(DoubleMatrix* A, DoubleMatrix* L, DoubleMatrix* U) {

    int n = (int)A->getNumRows();
    
    for (int j=0; j<n; j++)
        {
        for (int k=0; k<j; k++)
            for (int i=k+1; i<j; i++)
                (*A)(i,j) = (*A)(i,j) - (*A)(i,k) * (*A)(k,j);

        for (int k=0; k<j; k++)
            for (int i=j; i<n; i++)
                (*A)(i,j) = (*A)(i,j) - (*A)(i,k) * (*A)(k,j);

        for (int m=j+1; m<n; m++)
            (*A)(m,j) /= (*A)(j,j);
        }

    for (int row=0; row<n; row++)
        {
        for (int col=0; col<n; col++)
            {
            if ( row <= col )
                {
                (*U)(row,col) = (*A)(row,col);
                (*L)(row,col) = (row == col ? 1.0 : 0.0);
                }
            else
                {
                (*L)(row,col) = (*A)(row,col);
                (*U)(row,col) = 0.0;
                }
            }
        }
}

void MatrixMath::divideMatrixByPowerOfTwo(DoubleMatrix* M, int power) {

#   if defined(ROW_MAJOR_ORDER)

    int divisor = 1;
    for (int i=0; i<power; i++)
        divisor = divisor * 2;
    double factor = 1.0 / divisor;
    for (auto mPtr=M->begin(); mPtr != M->end(); mPtr++)
        *mPtr *= factor;

#   else

    int divisor = 1;
    for (int i=0; i<power; i++)
        divisor = divisor * 2;
        
    int nr = (int)M->getNumRows();
    int nc = (int)M->getNumCols();
        
    double factor = 1.0 / divisor;
    for (int i=0; i<nr; i++)
        for (int j=0; j<nc; j++)
            (*M)[i][j] *= factor;
            
#   endif
}

void MatrixMath::forwardSubstitutionRow(DoubleMatrix* L, double* b) {

    int n = (int)L->getNumRows();
    
    b[0] = b[0] / (*L)(0,0);
    for (int i=1; i<n; i++)
        {
        double dotProduct = 0.0;
        for (int j=0; j<i; j++)
            dotProduct += (*L)(i,j) * b[j];
        b[i] = (b[i] - dotProduct) / (*L)(i,i);
        }
}

void MatrixMath::gaussianElimination(DoubleMatrix* A, DoubleMatrix* B, DoubleMatrix* X, DoubleMatrix* L, DoubleMatrix* U, double* b) {

    int n = (int)A->getNumRows();

    computeLandU(A, L, U);

    for (int k=0; k<n; k++)
        {
        for (int i=0; i<n; i++)
            b[i] = (*B)(i,k);

        /* Answer of Ly = b (which is solving for y) is copied into b. */
        forwardSubstitutionRow(L, b);

        /* Answer of Ux = y (solving for x and the y was copied into b above)
           is also copied into b. */
        backSubstitutionRow(U, b);
        for (int i=0; i<n; i++)
            (*X)(i,k) = b[i];
        }
}

void MatrixMath::multiplicationByScalar(DoubleMatrix* M, double c) {

#   if defined(ROW_MAJOR_ORDER)
        
    for (auto mPtr=M->begin(); mPtr != M->end(); mPtr++)
        *mPtr *= c;
        
#   else

    int nr = (int)M->getNumRows();
    int nc = (int)M->getNumCols();
    for (int i=0; i<nr; i++)
        for (int j=0; j<nc; j++)
            (*M)[i][j] *= c;
            
#   endif
}

void MatrixMath::multiplicationByScalar(DoubleMatrix* M, double c, DoubleMatrix* Res) {

#   if defined(ROW_MAJOR_ORDER)
    
    if (M->size() != Res->size() || M->getNumRows() != Res->getNumRows())
        Msg::error("Error in matrix dimensions in multiplicationByScalar");
                
    for (auto mPtr=M->begin(), resPtr=Res->begin(); mPtr != M->end(); mPtr++, resPtr++)
        *resPtr = (*mPtr) * c;

#   else

    int nr = (int)M->getNumRows();
    int nc = (int)M->getNumCols();
    
    if (nr != Res->getNumRows() || nc != Res->getNumCols())
        Msg::error("Error in matrix dimensions in multiplicationByScalar");

    for (int i=0; i<nr; i++)
        for (int j=0; j<nc; j++)
            (*Res)[i][j] = (*M)[i][j] * c;
            
#   endif
}

void MatrixMath::multiplyTwoMatrices(DoubleMatrix* A, DoubleMatrix* B, DoubleMatrix* C, DoubleMatrix* temp) {

    if (A->getNumCols() != B->getNumRows())
        Msg::error("Problem with matrix dimensions in multiplyTwoMatrices");
        
    int nrA = (int)A->getNumRows();
    int ncA = (int)A->getNumCols();
    int ncB = (int)B->getNumCols();
    
    for (int i=0; i<nrA; i++)
        {
        for (int j=0; j<ncB; j++)
            {
            (*temp)(i,j) = 0.0;
            for (int k=0; k<ncA; k++)
                (*temp)(i,j) += (*A)(i,k) * (*B)(k,j);
            }
        }
    for (int i=0; i<nrA; i++)
        for (int j=0; j<ncB; j++)
            (*C)(i,j) = (*temp)(i,j);
}

void MatrixMath::setIdentity(DoubleMatrix* M) {

#   if defined(ROW_MAJOR_ORDER)

    int nr = (int)M->getNumRows();
    if (nr != M->getNumCols())
        Msg::error("Expectiung a square matrix in setIdentity");

    M->fill(0.0);
    for (int i=0; i<nr; i++)
        (*M)(i,i) = 1.0;

#   else

    int nr = (int)M->getNumRows();
    if (nr != M->getNumCols())
        Msg::error("Expectiung a square matrix in setIdentity");
        
    for (int i=0; i<nr; i++)
        for (int j=0; j<nr; j++)
            (*M)[i][j] = 0.0;
    for (int i=0; i<nr; i++)
        (*M)[i][i] = 1.0;
        
#   endif
}
