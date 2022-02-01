#ifndef Container_H
#define Container_H

#include <iostream>

template<typename T> class BufferTemplate {
    public:
               ~BufferTemplate();
        size_t size() { return Elements; }
        void   clear();
        bool   operator==(const BufferTemplate<T>& a) const;
        bool   operator!=(const BufferTemplate<T>& a) const;
        void   operator|=(T c);
        void   operator+=(T c);
        void   operator-=(T c);
        void   operator*=(T c);
        void   operator/=(T c);
        void   operator^=(T c);

    protected:
        BufferTemplate();
        BufferTemplate(const BufferTemplate<T> &a);
        void create(size_t elements);

        T*     Buffer;
        T*     EndBuffer;
        size_t Elements;
};

template<typename T> class ArrayTemplate: public BufferTemplate<T> {
    public:
             ArrayTemplate();
             explicit ArrayTemplate(size_t size);
             explicit ArrayTemplate(const ArrayTemplate<T>& a);
             void create(size_t elements) {__super::create(elements);}
        T    operator[](size_t i) const {return this.Buffer[i];}
        void operator+=(const ArrayTemplate<T>& a);
        void operator-=(const ArrayTemplate<T>& a);
};

template<typename T> class MatrixTemplate : public BufferTemplate<T> {
    public:
               MatrixTemplate();
               explicit MatrixTemplate(size_t rows, size_t cols);
               explicit MatrixTemplate(const MatrixTemplate<T>& a);
               void   create(size_t rows, size_t cols);
        size_t rows() { return Rows; }
        size_t cols() { return Cols; }
        T*     operator[](size_t i) const { return this.Buffer + i * Rows; }
        T      getValue(size_t i, size_t j) const { return this.Buffer[i * Rows + j]; }
        void   setValue(size_t i, size_t j, T value) {this.Buffer[i * Rows + j] = value; }
        void   transpose();
        bool   operator==(const MatrixTemplate<T>& a) const;
        bool   operator!=(const MatrixTemplate<T>& a) const;

    private:
        size_t Rows,
               Cols;
};

typedef ArrayTemplate<int>    IntArray;
typedef ArrayTemplate<double> DoubleArray;

typedef MatrixTemplate<int>    IntMatrix;
typedef MatrixTemplate<double> DoubleMatrix;

#endif