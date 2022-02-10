#ifndef Container_H
#define Container_H

#include <iostream>
#include <string>

#define ForElements(v) for (T* v = buffer; v < endBuffer; ++v)
#define ForLeftRight(a) const T* right = a; for (T* left = buffer; left < endBuffer; ++left, ++right)
#define IFVERIFY(condition) _ASSERT(condition); if (condition)



#pragma mark - BufferTemplate Definition -

template<typename T>
class BufferTemplate {

    public:
                   ~BufferTemplate(void);
        const T*    begin(void) const { return buffer; }
        T*          begin(void) { return buffer; }
        const T*    end(void) const { return endBuffer; }
        T*          end(void) { return endBuffer; }
        void        fill(T value);
        size_t      size(void) { return numElements; }
        void        setZero(void) { fill(0); }
        bool        operator==(const BufferTemplate<T>& a) const;
        bool        operator!=(const BufferTemplate<T>& a) const;
        void        operator+=(const BufferTemplate<T>& a);
        void        operator-=(const BufferTemplate<T>& a);
        void        operator|=(T c);
        void        operator+=(T c);
        void        operator-=(T c);
        void        operator*=(T c);
        void        operator/=(T c);
        void        operator^=(T c);

    protected:
                    BufferTemplate(void) = delete;
        explicit    BufferTemplate(const BufferTemplate<T> &a);
        explicit    BufferTemplate(size_t ne);
        void        create(size_t ne);
        T*          buffer;
        T*          endBuffer;
        size_t      numElements;
};



#pragma mark - ArrayTemplate Definition -

template<typename T>
class ArrayTemplate: public BufferTemplate<T> {

    public:
                            ArrayTemplate(void);
        explicit            ArrayTemplate(size_t size);
        explicit            ArrayTemplate(const ArrayTemplate<T>& a);
        ArrayTemplate<T>&   operator=(const ArrayTemplate<T>& rhs);
        void                create(size_t elements) { BufferTemplate<T>::create(elements); }
        T                   operator[](size_t i) const { return this->buffer[i]; }
        T                   getValue(size_t i) const { return this->buffer[i]; }
        
    private:
        void                copy(const ArrayTemplate<T>& other);
};



#pragma mark - MatrixTemplate Definition -

template<typename T>
class MatrixTemplate : public BufferTemplate<T> {

    public:
                            MatrixTemplate(void);
        explicit            MatrixTemplate(size_t rows, size_t cols);
        explicit            MatrixTemplate(const MatrixTemplate<T>& m);
        MatrixTemplate<T>&  operator=(const MatrixTemplate<T>& rhs);
        T&                  operator()(size_t r, size_t c) { return this->buffer[r * numCols + c]; }
        const T&            operator()(size_t r, size_t c) const { return this->buffer[r * numCols + c]; }
        bool                operator==(const MatrixTemplate<T>& m) const;
        bool                operator!=(const MatrixTemplate<T>& m) const;
        void                create(size_t rows, size_t cols);
        size_t              rows(void) const { return numRows; }
        size_t              cols(void) const { return numCols; }
        T                   getValue(size_t r, size_t c) const;
        void                setValue(size_t r, size_t c, T value);
        void                setIdentity(T value=1);
        void                print(void);
        void                print(std::string title);
        void                multiply(const MatrixTemplate<T>& m, MatrixTemplate<T>& result) const;
        void                transpose(MatrixTemplate<T>& result);

    protected:
        void                copy(const MatrixTemplate<T>& other);
        size_t              numRows,
                            numCols;
};



#pragma mark - BufferTemplate Implementation -

template<typename T>
BufferTemplate<T>::BufferTemplate(const BufferTemplate<T>& a) {

    create(a.numElements);
    if (numElements)
        memcpy(buffer, a.buffer, numElements * sizeof(T));
}

template<typename T>
BufferTemplate<T>::BufferTemplate(size_t ne) {

    numElements = 0;
    create(ne);
}

template<typename T>
BufferTemplate<T>::~BufferTemplate(void) {

    delete [] buffer;
}

template<typename T>
void BufferTemplate<T>::create(size_t ne) {

    if (ne != numElements)
        {
        if (numElements > 0)
            delete [] buffer;
        if (ne > 0)
            buffer = new T[ne];
        else
            buffer = NULL;
        }
    numElements = ne;
    endBuffer = buffer + numElements;
}

template<typename T>
void BufferTemplate<T>::fill(T value) {

    memset(buffer, value, numElements * sizeof(T));
}

template<typename T>
void BufferTemplate<T>::operator+=(T c) {

    ForElements(e)
        *e += c;
}

template<typename T>
void BufferTemplate<T>::operator-=(T c) {

    ForElements(e)
        *e -= c;
}

template<typename T>
void BufferTemplate<T>::operator*=(T c) {

    ForElements(e)
        *e *= c;
}

template<typename T>
void BufferTemplate<T>::operator/=(T c) {

    ForElements(e)
        *e /= c;
}

template<typename T>
void BufferTemplate<T>::operator|=(T c) {

    ForElements(e)
        *e |= c;
}

template<typename T>
void BufferTemplate<T>::operator^=(T c) {

    ForElements(e)
        *e ^= c;
}

template<typename T>
bool BufferTemplate<T>::operator==(const BufferTemplate<T>& a) const {

    if (numElements != a.numElements)
        return false;
    ForLeftRight(a)
        {
        if (*left != *right)
            return false;
        }
    return true;
}

template<typename T>
bool BufferTemplate<T>::operator!=(const BufferTemplate<T>& a) const {

    if (numElements != a.numElements)
        return true;
     ForLeftRight(a)
        {
        if (*left != *right)
            return true;
        }
    return false;
}

template<typename T>
void BufferTemplate<T>::operator+=(const BufferTemplate<T>& a) {

    ForLeftRight(a)
        *left += *right;
}

template<typename T>
void BufferTemplate<T>::operator-=(const BufferTemplate<T>& a) {

    ForLeftRight(a)
        *left -= *right;
}



#pragma mark - ArrayTemplate Implementation -

template<typename T>
ArrayTemplate<T>::ArrayTemplate(void) {

    create(1);
}

template<typename T>
ArrayTemplate<T>::ArrayTemplate(size_t ne) {

    create(ne);
}

template<typename T>
ArrayTemplate<T>::ArrayTemplate(const ArrayTemplate<T>& a) {

    copy(a);
}

template<typename T>
ArrayTemplate<T>& ArrayTemplate<T>::operator=(const ArrayTemplate<T>& rhs) {

    if (this != &rhs)
        copy(rhs);
    return *this;
}

template<typename T>
void ArrayTemplate<T>::copy(const ArrayTemplate<T>& other) {

    create(other.numElements);
    memcpy(this->buffer, other.buffer, sizeof(T) * this->numElements);
}


#pragma mark - MatrixTemplate Implementation -

template<typename T>
MatrixTemplate<T>::MatrixTemplate(void) : BufferTemplate<T>::BufferTemplate(1) {

    numRows = 1;
    numCols = 1;
}

template<typename T>
MatrixTemplate<T>::MatrixTemplate(const MatrixTemplate<T>& m) : BufferTemplate<T>::BufferTemplate(m) {

    copy(m);
}

template<typename T>
MatrixTemplate<T>::MatrixTemplate(size_t nr, size_t nc) : BufferTemplate<T>::BufferTemplate(nr * nc) {

    numRows = nr;
    numCols = nc;
    BufferTemplate<T>::fill(0);
}

template<typename T>
MatrixTemplate<T>& MatrixTemplate<T>::operator=(const MatrixTemplate<T>& rhs) {

    if (this != &rhs)
        copy(rhs);
    return *this;
}

template<typename T>
void MatrixTemplate<T>::print(void) {

    std::cout << "{";
    for (int i=0; i<numRows; i++)
        {
        std::cout << "{";
        for (int j=0; j<numCols; j++)
            {
            std::cout << getValue(i, j);
            if (j != numCols-1)
                std::cout << ",";
            }
        if (i != numRows-1)
            std::cout << "}," << std::endl;
        }
    std::cout << "}}" << std::endl;
}

template<typename T>
void MatrixTemplate<T>::print(std::string title) {

    std::cout << title << std::endl;
    print();
}

template<typename T>
void MatrixTemplate<T>::copy(const MatrixTemplate<T>& other) {

    create(other.numRows, other.numCols);
    numRows = other.numRows;
    numCols = other.numCols;
    memcpy(this->buffer, other.buffer, numRows * numCols * sizeof(T));
}

template<typename T>
void MatrixTemplate<T>::create(size_t nr, size_t nc) {

    BufferTemplate<T>::create(nr * nc);
    numRows = nr;
    numCols = nc;
}

template<typename T> T
MatrixTemplate<T>::getValue(size_t r, size_t c) const {

    return this->buffer[r * numCols + c];
}

template<typename T>
void MatrixTemplate<T>::setValue(size_t r, size_t c, T value) {

    this->buffer[r * numCols + c] = value;
}

template<typename T>
void MatrixTemplate<T>::setIdentity(T value) {

    IFVERIFY(numRows > 0 && numRows == numCols)
        {
        BufferTemplate<T>::fill(0);
        for (int i = 0; i < numRows; ++i)
            setValue(i, i, value);
        }
}

template<typename T>
bool MatrixTemplate<T>::operator==(const MatrixTemplate<T>& m) const {

    if (numRows != m.numRows || numCols != m.numCols)
        return false;
    return BufferTemplate<T>::operator==(m);
}

template<typename T>
bool MatrixTemplate<T>::operator!=(const MatrixTemplate<T>& m) const {

    if (numRows != m.numRows || numCols != m.numCols)
        return true;
    return BufferTemplate<T>::operator!=(m);
}

template<typename T>
void MatrixTemplate<T>::transpose(MatrixTemplate<T>& result)  {

    result.Create(numCols, numRows);
    for (int r = 0; r < numRows; ++r)
        {
        for (int c = 0; c < numCols; ++c)
            result.setValue(c, r, getValue(r, c));
        }
}

template<typename T>
void MatrixTemplate<T>::multiply(const MatrixTemplate<T>& m, MatrixTemplate<T>& result) const {

    IFVERIFY (numCols == m.numRows)
        {
        result.create(numRows, m.numCols);
        for (size_t i = 0; i < numRows; ++i)
            {
            for (size_t j = 0; j < m.numCols; ++j)
                {
                T total = getValue(i, 0) * m.getValue(0, j);
                for (size_t z = 1; z < m.numRows; ++z)
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
