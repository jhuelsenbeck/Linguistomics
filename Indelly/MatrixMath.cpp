#include "MatrixMath.hpp"
#include "Msg.hpp"

MathCache::MathCache(int numStates):
    scratch1(numStates, numStates),
    scratch2(numStates, numStates),
    scratchmult(numStates, numStates)
{
    scratchVec = new double[numStates];
}

MathCache::~MathCache() {
    delete scratchVec;
}

void MathCache::backSubstitutionRow(DoubleMatrix& U, double* b) {

    int n = (int)U.getNumRows();
    int np1 = n + 1;

    auto bp = b + n - 1;

    auto urow  = U.end() - n;
    auto u     = U.end() - 1;
    auto udiag = u;

    *bp-- /= *u;


    for (int i = n - 2; i >= 0; i--)
    {
        u     -= n;
        urow  -= n;
        udiag -= np1;

        int j = i + 1;

        double dotProduct = 0.0;
        for (auto uj = urow + j, end = urow + n, bj = b + j; uj < end; uj++, bj++)
            dotProduct += *uj * *bj;

        *bp = (*bp - dotProduct) / *udiag;
        --bp;
    }
}

void MathCache::forwardSubstitutionRow(DoubleMatrix& L, double* const b) {
    auto lrows = L.getNumRows();
    auto lcols = L.getNumCols();
    auto lcols1 = lcols+1;
    auto lrow = L.begin();
    auto ldiag = lrow;
    auto bp = b;

    *bp++ /= *lrow;
    for (size_t i = 1; i < lrows; i++)
    {
        lrow += lcols;
        ldiag += lcols1;

        double dotProduct = 0.0;
        for (auto lp = lrow, lend = lrow + i, bj = b; lp < lend; ++lp, ++bj)
            dotProduct += *lp * *bj;

        *bp = (*bp - dotProduct) / *ldiag;
        ++bp;
    }
}

void MathCache::computeLandU(DoubleMatrix& A, DoubleMatrix& L, DoubleMatrix& U) {

    size_t n = A.getNumRows();

    for (size_t j = 0; j < n; j++)
    {
        for (size_t k = 0; k < j; k++)
            for (size_t i = k + 1; i < j; i++)
                A(i, j) -= A(i, k) * A(k, j);

        for (size_t k = 0; k < j; k++)
            for (size_t i = j; i < n; i++)
                A(i, j) -= A(i, k) * A(k, j);

        for (size_t m = j + 1; m < n; m++)
            A(m, j) /= A(j, j);
    }

    auto u = U.begin();
    auto l = L.begin();
    auto a = A.begin();
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

void MathCache::gaussianElimination(DoubleMatrix& A, DoubleMatrix& B, DoubleMatrix& X, DoubleMatrix& L) {

    auto n  = A.getNumRows();
    auto b  = scratchVec;

    computeLandU(A, L, scratch2);

    auto brow = B.begin();
    auto xrow = X.begin();

    for (int k = 0; k < n; k++)
    {
        for (auto bp = b, bend = b + n, br = brow; bp < bend; ++bp, br += n)
            *bp = *br;

        /* Answer of Ly = b (which is solving for y) is copied into b. */
        forwardSubstitutionRow(L, b);

        /* Answer of Ux = y (solving for x and the y was copied into b above)
           is also copied into b. */
        backSubstitutionRow(scratch2, b);

        for (auto bp = b, bend = b + n, xr = xrow; bp < bend; ++bp, xr += n)
            *xr = *bp;

        ++brow;
        ++xrow;
    }
}

void MathCache::multiply(DoubleMatrix& A, DoubleMatrix&B) {
    A.multiply(B, scratchmult);
    B.copy(scratchmult);
}
