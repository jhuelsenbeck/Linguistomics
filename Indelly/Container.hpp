#ifndef Container_H
#define Container_H

#include <iostream>
#define ForElements(v) for (T* v = Buffer; v < endBuffer; ++v)
#define ForLeftRight(a) const T* right = a; for (T* left = Buffer; left < endBuffer; ++left, ++right)
#define IFVERIFY(condition) _ASSERT(condition); if (condition)


#pragma mark - BufferTemplate Definition -

template<typename T> class BufferTemplate {

    public:
               ~BufferTemplate(void);
        size_t  size(void) { return numElements; }
        void    setZero(void) { fill(0); }
        void    fill(T value);
        bool    operator==(const BufferTemplate<T>& a) const;
        bool    operator!=(const BufferTemplate<T>& a) const;
        void    operator+=(const BufferTemplate<T>& a);
        void    operator-=(const BufferTemplate<T>& a);
        void    operator|=(T c);
        void    operator+=(T c);
        void    operator-=(T c);
        void    operator*=(T c);
        void    operator/=(T c);
        void    operator^=(T c);

    protected:
                BufferTemplate(void);
                BufferTemplate(const BufferTemplate<T> &a);
                BufferTemplate(size_t elements);
        void    create(size_t elements);
        T*      Buffer;
        T*      endBuffer;
        size_t  numElements;
};

#pragma mark - ArrayTemplate Definition -

template<typename T> class ArrayTemplate: public BufferTemplate<T> {

    public:
                ArrayTemplate(void);
                explicit ArrayTemplate(size_t size);
                explicit ArrayTemplate(const ArrayTemplate<T>& a); // copy
                void create(size_t elements) { BufferTemplate<T>::create(elements);}
        T       operator[](size_t i) const {return this->Buffer[i];}
        T       getValue(size_t i) const { return this->Buffer[i]; }
};

#pragma mark - MatrixTemplate Definition -

template<typename T> class MatrixTemplate : public BufferTemplate<T> {

    public:
                MatrixTemplate(void);
                explicit MatrixTemplate(size_t rows, size_t cols);
                explicit MatrixTemplate(const MatrixTemplate<T>& m);  // copy
        void    create(size_t rows, size_t cols);
        size_t  rows(void) const { return Rows; }
        size_t  cols(void) const { return Cols; }
        T       getValue(size_t r, size_t c) const;
        bool    isSquare(void) { return (Rows == Cols); } 
        void    setValue(size_t r, size_t c, T value);
        void    setIdentity(T value=1);
        bool    operator==(const MatrixTemplate<T>& m) const;
        bool    operator!=(const MatrixTemplate<T>& m) const;
        void    multiply(const MatrixTemplate<T>& m, MatrixTemplate<T>& result) const;
        void    transpose(MatrixTemplate<T>& result);

    private:
        size_t  Rows,
                Cols;
};



#pragma mark - BufferTemplate Implementation -

template<typename T> BufferTemplate<T>::BufferTemplate(void) {

    numElements = 0;
    Buffer = NULL;
    endBuffer = NULL;
}

template<typename T> BufferTemplate<T>::BufferTemplate(const BufferTemplate<T>& a) {

    create(a.numElements);
    if (numElements)
        memcpy(Buffer, a.Buffer, numElements * sizeof(T));
}

template<typename T> BufferTemplate<T>::BufferTemplate(size_t elements) {

    create(elements);
}

template<typename T> BufferTemplate<T>::~BufferTemplate(void) {

    delete [] Buffer;
}

template<typename T> void BufferTemplate<T>::create(size_t ne) {

    if (ne > numElements)
        {
        delete[] Buffer;
        if (ne > 0)
            Buffer = new T[ne];
        else
            Buffer = NULL;
        }
    numElements = ne;
    endBuffer = Buffer + numElements;
}

template<typename T> void BufferTemplate<T>::fill(T value) {

    memset(Buffer, value, numElements);
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

    if (numElements != a.numElements)
        return false;
    ForLeftRight(a)
        {
        if (*left != *right)
            return false;
        }
    return true;
}

template<typename T> bool BufferTemplate<T>::operator!=(const BufferTemplate<T>& a) const {

    if (numElements != a.numElements)
        return true;
     ForLeftRight(a)
        {
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

#pragma mark - MatrixTemplate Implementation -

template<typename T> MatrixTemplate<T>::MatrixTemplate(void) {

    Rows = 0;
    Cols = 0;
}

template<typename T> MatrixTemplate<T>::MatrixTemplate(const MatrixTemplate<T>& m) : BufferTemplate<T>::BufferTemplate(m) {

    Rows = m.Rows;
    Cols = m.Cols;
}

template<typename T>
MatrixTemplate<T>::MatrixTemplate(size_t rows, size_t cols) : BufferTemplate<T>::BufferTemplate(rows * cols) {

    Rows = rows;
    Cols = cols;
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

    IFVERIFY (Cols == m.Rows)
        {
        result.create(Rows, m.Cols);
        for (size_t i = 0; i < Rows; ++i)
            {
            for (size_t j = 0; j < m.Cols; ++j)
                {
                T total = getValue(i, 0) * m.getValue(0, j);
                for (size_t z = 1; z < m.Rows; ++z)
                    total += getValue(i, z) * m.getValue(z, j);
                result.SetValue(i, j, total);
                }
            }
        }
}

#pragma mark - Type Definitions -

typedef ArrayTemplate<int>    IntArray;
typedef ArrayTemplate<double> DoubleArray;

typedef MatrixTemplate<int>    IntMatrix;
typedef MatrixTemplate<double> MyDoubleMatrix;

#endif
