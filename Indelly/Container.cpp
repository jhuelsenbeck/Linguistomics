#include "Container.hpp"

//==========================================================================

#define ForElements(v) for (T* v = this.Buffer; v < this.EndBuffer; ++v)
#define ForLeftRight(a) const T* right = a; for (T* left = this.Buffer; left < this.EndBuffer; ++left, ++right)


template<typename T> BufferTemplate<T>::BufferTemplate() {
    Elements = 0;
    Buffer = NULL;
    EndBuffer = NULL;
}

template<typename T> BufferTemplate<T>::~BufferTemplate() {
    delete[] Buffer;
}

template<typename T> void BufferTemplate<T>::create(size_t elements) {
    if (elements > Elements)
        {
        delete[] Buffer;
        Buffer = new T[elements];
        }
    Elements = elements;
    EndBuffer = Buffer + elements;
}


template<typename T> void BufferTemplate<T>::clear() {
    memset(Buffer, 0, Elements); 
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
        if (*left == *right++)
            return false;
    }
    return true;
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
//==========================================================================

template<typename T> ArrayTemplate<T>::ArrayTemplate() {
}

template<typename T> ArrayTemplate<T>::ArrayTemplate(size_t size) {
    create(size);
}

template<typename T> void ArrayTemplate<T>::operator+=(const ArrayTemplate<T>& a) {
    ForLeftRight(a)
        *left += *right++;
}

template<typename T> void ArrayTemplate<T>::operator-=(const ArrayTemplate<T>& a) {
    ForLeftRight(a)
        *left -= *right++;
}
//==========================================================================

template<typename T> MatrixTemplate<T>::MatrixTemplate() {
    Rows = 0;
    Cols = 0;
}

template<typename T> MatrixTemplate<T>::MatrixTemplate(size_t rows, size_t cols) {
    create(rows, cols);
}

template<typename T> void MatrixTemplate<T>::create(size_t rows, size_t cols) {
    Rows = rows;
    Cols = cols;
    __super::create(rows * cols);
}

template<typename T> bool MatrixTemplate<T>::operator==(const MatrixTemplate<T>& a) const {
    if (Rows != a.Rows || Cols != a.Cols)
        return false;
    return __super::operator==(a);
}

template<typename T> bool MatrixTemplate<T>::operator!=(const MatrixTemplate<T>& a) const {
    if (Rows != a.Rows || Cols != a.Cols)
        return true;
    return __super::operator!=(a);
}
//==========================================================================
