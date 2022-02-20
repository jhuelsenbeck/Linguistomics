#include "MatrixMath.hpp"
#include "Msg.hpp"

MathCache::MathCache(int numStates):
    scratch1(numStates, numStates),
    scratch2(numStates, numStates)
{
    scratchVec = new double[numStates];
}

MathCache::~MathCache() {
    delete scratchVec;
}




void MatrixMath::backSubstitutionRow(DoubleMatrix* U, double* b) {

    int n = (int)U->getNumRows();

    b[n - 1] /= (*U)(n - 1, n - 1);
    for (int i = n - 2; i >= 0; i--)
    {
        double dotProduct = 0.0;
        for (int j = i + 1; j < n; j++)
            dotProduct += (*U)(i, j) * b[j];
        b[i] = (b[i] - dotProduct) / (*U)(i, i);
    }
}

void MatrixMath::forwardSubstitutionRow(DoubleMatrix* L, double* b) {

    size_t n = L->getNumRows();

    b[0] = b[0] / (*L)(0, 0);
    for (size_t i = 1; i < n; i++)
    {
        double dotProduct = 0.0;
        for (size_t j = 0; j < i; j++)
            dotProduct += (*L)(i, j) * b[j];
        b[i] = (b[i] - dotProduct) / (*L)(i, i);
    }
}

void MatrixMath::computeLandU(DoubleMatrix* A, DoubleMatrix* L, DoubleMatrix* U) {

    size_t n = A->getNumRows();

    for (size_t j = 0; j < n; j++)
    {
        for (size_t k = 0; k < j; k++)
            for (size_t i = k + 1; i < j; i++)
                (*A)(i, j) -= (*A)(i, k) * (*A)(k, j);

        for (size_t k = 0; k < j; k++)
            for (size_t i = j; i < n; i++)
                (*A)(i, j) -= (*A)(i, k) * (*A)(k, j);

        for (size_t m = j + 1; m < n; m++)
            (*A)(m, j) /= (*A)(j, j);
    }

    auto u = U->begin();
    auto l = L->begin();
    auto a = A->begin();
    for (size_t row = 0; row < n; row++)
    {
        for (size_t col = 0; col < n; col++)
        {
            if (row <= col)
            {
                *u = *a;
                *l = (row == col ? 1.0 : 0.0);
            }
            else
            {
                *l = *a;
                *u = 0.0;
            }

            ++l;
            ++u;
            ++a;
        }
    }
}

void MatrixMath::gaussianElimination(DoubleMatrix* A, DoubleMatrix* B, DoubleMatrix* X, DoubleMatrix* L, DoubleMatrix* U, double* b) {

    auto n = A->getNumRows();

    computeLandU(A, L, U);

    for (int k = 0; k < n; k++)
    {
        for (int i = 0; i < n; i++)
            b[i] = (*B)(i, k);

        /* Answer of Ly = b (which is solving for y) is copied into b. */
        forwardSubstitutionRow(L, b);

        /* Answer of Ux = y (solving for x and the y was copied into b above)
           is also copied into b. */
        backSubstitutionRow(U, b);
        for (int i = 0; i < n; i++)
            (*X)(i, k) = b[i];
    }
}
