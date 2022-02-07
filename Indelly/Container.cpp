#include "Container.hpp"

//==========================================================================

#define ForElements(v) for (T* v = Buffer; v < EndBuffer; ++v)
#define ForLeftRight(a) const T* right = a; for (T* left = Buffer; left < EndBuffer; ++left, ++right)

// not sure what the platform neutral version of this is
#define IFVERIFY(condition) _ASSERT(condition); if (condition)

// John, use this macro below if the above macro does not compile
//#define IFVERIFY(condition) if (condition)


template<typename T> BufferTemplate<T>::BufferTemplate() {
    Elements = 0;
    Buffer = NULL;
    EndBuffer = NULL;
}

template<typename T> BufferTemplate<T>::BufferTemplate(const BufferTemplate<T>& a) {
    create(a.Elements);
    if (Elements)
        memcpy(Buffer, a.Buffer, Elements * sizeof(T));
}

template<typename T> BufferTemplate<T>::BufferTemplate(size_t elements) {
    create(elements);
}

template<typename T> BufferTemplate<T>::~BufferTemplate() {
    delete[] Buffer;
}

template<typename T> void BufferTemplate<T>::create(size_t elements) {
    if (elements > Elements)
        {
        delete[] Buffer;
        if (elements > 0)
            Buffer = new T[elements];
        else
            Buffer = NULL;
        }
    Elements = elements;
    EndBuffer = Buffer + elements;
}

template<typename T> void BufferTemplate<T>::fill(T value) {
    memset(Buffer, value, Elements); 
}

template<typename T> void BufferTemplate<T>::operator+=(T c) {
    ForElements(e)
        *e += c;
}

template<typename T> void BufferTemplate<T>::operator-=(T c) {
    ForElements(e)
        *e -= c;
}

template<typename T> void BufferTemplate<T>::operator*=(T c) {
    ForElements(e)
        *e *= c;
}

template<typename T> void BufferTemplate<T>::operator/=(T c) {
    ForElements(e)
        *e /= c;
}

template<typename T> void BufferTemplate<T>::operator|=(T c) {
    ForElements(e)
        *e |= c;
}

template<typename T> void BufferTemplate<T>::operator^=(T c) {
    ForElements(e)
        *e ^= c;
}

template<typename T> bool BufferTemplate<T>::operator==(const BufferTemplate<T>& a) const {
    if (Elements != a.Elements) 
        return false;
    ForLeftRight(a) 
        {
        if (*left != *right)
            return false;
        }
    return true;
}

template<typename T> bool BufferTemplate<T>::operator!=(const BufferTemplate<T>& a) const {
    if (Elements != a.Elements)
        return true;
     ForLeftRight(a) {
        if (*left != *right)
            return true;
    }
    return false;
}

template<typename T> void BufferTemplate<T>::operator+=(const BufferTemplate<T>& a) {
    ForLeftRight(a)
        *left += *right;
}

template<typename T> void BufferTemplate<T>::operator-=(const BufferTemplate<T>& a) {
    ForLeftRight(a)
        *left -= *right;
}

//==========================================================================

template<typename T> ArrayTemplate<T>::ArrayTemplate() {
}

template<typename T> ArrayTemplate<T>::ArrayTemplate(size_t size) {
    create(size);
}

template<typename T> ArrayTemplate<T>::ArrayTemplate(const ArrayTemplate<T>& a) {
    copy(a);
}
//==========================================================================

#define ForRows for (T* row = Buffer; row < EndBuffer; row += Cols)
#define ForColumns for (T* col = row; col < row + Cols; ++col)

template<typename T> MatrixTemplate<T>::MatrixTemplate() {
    Rows = 0;
    Cols = 0;
}

template<typename T>
MatrixTemplate<T>::MatrixTemplate(size_t rows, size_t cols):
    BufferTemplate<T>::BufferTemplate(rows * cols)
{
    Rows = rows;
    Cols = cols;
}

template<typename T> MatrixTemplate<T>::MatrixTemplate(const MatrixTemplate<T>& m):
    BufferTemplate<T>::BufferTemplate(m)
{
    Rows = m.Rows;
    Cols = m.Cols;
}

template<typename T> void MatrixTemplate<T>::create(size_t rows, size_t cols) {
    BufferTemplate<T>::create(rows * cols);
    Rows = rows;
    Cols = cols;
}

template<typename T> T MatrixTemplate<T>::getValue(size_t r, size_t c) const { 
    return this->Buffer[r * Cols + c];
}

template<typename T> void MatrixTemplate<T>::setValue(size_t r, size_t c, T value) { 
    this->Buffer[r * Cols + c] = value;
}

template<typename T> void MatrixTemplate<T>::setIdentity(T value) {
    IFVERIFY(Rows > 0 && Rows == Cols)
    {
        BufferTemplate<T>::fill(0);
        for (int i = 0; i < Rows; ++i)
            setValue(i, i, value);
    }
}

template<typename T> bool MatrixTemplate<T>::operator==(const MatrixTemplate<T>& m) const {
    if (Rows != m.Rows || Cols != m.Cols)
        return false;
    return BufferTemplate<T>::operator==(m);
}

template<typename T> bool MatrixTemplate<T>::operator!=(const MatrixTemplate<T>& m) const {
    if (Rows != m.Rows || Cols != m.Cols)
        return true;
    return BufferTemplate<T>::operator!=(m);
}

template<typename T> void MatrixTemplate<T>::transpose(MatrixTemplate<T>& result)  {
    result.Create(Cols, Rows);
    for (int r = 0; r < Rows; ++r) 
        {
        for (int c = 0; c < Cols; ++c)
            result.setValue(c, r, getValue(r, c));
        }
}

template<typename T> void MatrixTemplate<T>::multiply(const MatrixTemplate<T>& m, MatrixTemplate<T>& result) const {
    IFVERIFY (Cols == m.Rows) {
        result.create(Rows, m.Cols);
        for (size_t i = 0; i < Rows; ++i) {
            for (size_t j = 0; j < m.Cols; ++j) {
                T total = getValue(i, 0) * m.getValue(0, j);
                for (size_t z = 1; z < m.Rows; ++z)
                    total += getValue(i, z) * m.getValue(z, j);
                result.SetValue(i, j, total);
            }
        }
    }
}

//==========================================================================
