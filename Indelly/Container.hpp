#ifndef Container_H
#define Container_H

#include <iostream>

template<typename T> class BufferTemplate {
    public:
               ~BufferTemplate();
        size_t size() { return Elements; }
        void   setZero() { fill(0);}
        void   fill(T value);
        bool   operator==(const BufferTemplate<T>& a) const;
        bool   operator!=(const BufferTemplate<T>& a) const;
        void   operator+=(const BufferTemplate<T>& a);
        void   operator-=(const BufferTemplate<T>& a);
        void   operator|=(T c);
        void   operator+=(T c);
        void   operator-=(T c);
        void   operator*=(T c);
        void   operator/=(T c);
        void   operator^=(T c);

    protected:
             BufferTemplate();
             BufferTemplate(const BufferTemplate<T> &a);  // copy
             BufferTemplate(size_t elements);
        void create(size_t elements);

        T*     Buffer;
        T*     EndBuffer;
        size_t Elements;
};

template<typename T> class ArrayTemplate: public BufferTemplate<T> {
    public:
             ArrayTemplate();
             explicit ArrayTemplate(size_t size);
             explicit ArrayTemplate(const ArrayTemplate<T>& a); // copy
             void create(size_t elements) { BufferTemplate<T>::create(elements);}
        T    operator[](size_t i) const {return this->Buffer[i];}
        T    getValue(size_t i) const { return this->Buffer[i]; }
};

template<typename T> class MatrixTemplate : public BufferTemplate<T> {
    public:
               MatrixTemplate();
               explicit MatrixTemplate(size_t rows, size_t cols);
               explicit MatrixTemplate(const MatrixTemplate<T>& m);  // copy
        void   create(size_t rows, size_t cols);
        size_t rows() const { return Rows; }
        size_t cols() const { return Cols; }
        T      getValue(size_t r, size_t c) const;
        void   setValue(size_t r, size_t c, T value);
        void   setIdentity(T value=1);
        bool   operator==(const MatrixTemplate<T>& m) const;
        bool   operator!=(const MatrixTemplate<T>& m) const;
        void   multiply(const MatrixTemplate<T>& m, MatrixTemplate<T>& result) const;
        void   transpose(MatrixTemplate<T>& result);

    private:
        size_t Rows,
               Cols;
};

typedef ArrayTemplate<int>    IntArray;
typedef ArrayTemplate<double> DoubleArray;

typedef MatrixTemplate<int>    IntMatrix;
typedef MatrixTemplate<double> DoubleMatrix;

#endif
