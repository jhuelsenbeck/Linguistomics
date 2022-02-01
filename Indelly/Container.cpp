#include "Container.hpp"

//==========================================================================

#define ForElements(v) for (T* v = this.Buffer; v < this.EndBuffer; ++v)
#define ForLeftRight(a) const T* right = a; for (T* left = this.Buffer; left < this.EndBuffer; ++left, ++right)


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


template<typename T> void BufferTemplate<T>::clear() {
    memset(Buffer, 0, Elements); 
}

template<typename T> void BufferTemplate<T>::operator+=(T c) {
    ForElements(e)
        * e += c;
}

template<typename T> void BufferTemplate<T>::operator-=(T c) {
    ForElements(e)
        * e -= c;
}

template<typename T> void BufferTemplate<T>::operator*=(T c) {
    ForElements(e)
        * e *= c;
}

template<typename T> void BufferTemplate<T>::operator/=(T c) {
    ForElements(e)
        * e /= c;
}

template<typename T> void BufferTemplate<T>::operator|=(T c) {
    ForElements(e)
        * e |= c;
}

template<typename T> void BufferTemplate<T>::operator^=(T c) {
    ForElements(e)
        * e ^= c;
}

template<typename T> bool BufferTemplate<T>::operator==(const BufferTemplate<T>& a) const {
    if (Elements != a.Elements) 
        return false;
    ForLeftRight(a) 
        {
        if (*left != *right++)
            return false;
        }
    return true;
}

template<typename T> bool BufferTemplate<T>::operator!=(const BufferTemplate<T>& a) const {
    if (Elements != a.Elements)
        return true;
     ForLeftRight(a) {
        if (*left != *right++)
            return true;
    }
    return false;
}

template<typename T> void BufferTemplate<T>::operator+=(const BufferTemplate<T>& a) {
    ForLeftRight(a)
        * left += *right++;
}

template<typename T> void BufferTemplate<T>::operator-=(const BufferTemplate<T>& a) {
    ForLeftRight(a)
        * left -= *right++;
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

template<typename T> MatrixTemplate<T>::MatrixTemplate() {
    Rows = 0;
    Cols = 0;
}

template<typename T> MatrixTemplate<T>::MatrixTemplate(size_t rows, size_t cols) {
    create(rows, cols);
}

template<typename T> MatrixTemplate<T>::MatrixTemplate(const MatrixTemplate<T>& m):
    __super::BufferTemplate<T>(m)
{
    Rows = m.Rows;
    Cols = m.Cols;
}

template<typename T> void MatrixTemplate<T>::create(size_t rows, size_t cols) {
    __super::create(rows * cols);
    Rows = rows;
    Cols = cols;
}

template<typename T> bool MatrixTemplate<T>::operator==(const MatrixTemplate<T>& m) const {
    if (Rows != m.Rows || Cols != m.Cols)
        return false;
    return __super::operator==(m);
}

template<typename T> bool MatrixTemplate<T>::operator!=(const MatrixTemplate<T>& m) const {
    if (Rows != m.Rows || Cols != m.Cols)
        return true;
    return __super::operator!=(m);
}

template<typename T> void MatrixTemplate<T>::transpose()  {
    if (Rows == 1 && Cols == 1)
      return;

    for (int r = 0; r < Rows; ++r) 
        {
        for (int c = 0; c < Cols; ++c)
            setValue(c, r, getValue(r, c));
        }
}

template<typename T> MatrixTemplate<T>& MatrixTemplate<T>::operator*(const MatrixTemplate<T>& m) const {
    if (Cols == m.Rows) {
        auto result = new MatrixTemplate<T>(Rows, m.Cols);
        for (size_t i = 0; i < Rows; ++i) {
            for (size_t j = 0; j < m.Cols; ++j) {
                T total = getValue(i, 0) * m.getValue(0, j);
                for (size_t z = 1; z < m.Rows; ++z)
                    total += getValue(i, z) * m.getValue(z, j);
                result.SetValue(i, j, total);
            }
        }
        return result;
    }
    _ASSERT("invalid multiplication");
    return *this;
}

//==========================================================================
