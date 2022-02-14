#ifndef Container_H
#define Container_H

#include <iomanip>
#include <iostream>
#include <string>

#define ForElements(v) for (auto v = this->begin(), end = this->end(); v < end; ++v)
#define ForLeftRight(a) const auto* right = a.begin(); for (auto left = this->begin(), end = this->end(); left < end; ++left, ++right)
#define SIZEASSERT(condition)

#pragma mark - BufferTemplate Definition -


template<typename T>
class BufferTemplate {

    public:
        explicit BufferTemplate(const BufferTemplate<T>& a) {
            create(a.numElements);
            if (numElements)
                memcpy(buffer, a.buffer, numElements * sizeof(T));
        }

        explicit BufferTemplate(size_t ne) {

            numElements = 0;
            create(ne);
        }

        ~BufferTemplate(void) {

            delete[] buffer;
        }

        
        void create(size_t ne) {

            if (ne != numElements)
            {
                if (numElements > 0)
                    delete[] buffer;
                if (ne > 0)
                    buffer = new T[ne];
                else
                    buffer = NULL;
            }
            numElements = ne;
            endBuffer = buffer + numElements;
        }

        void setZero() {
            memset(this->begin(), 0, numElements * sizeof(T));
        }

        void copy(const BufferTemplate<T> &a) {

            create(a.numElements);
            memcpy(this->begin(), a.begin(), sizeof(T)* this->numElements);
        }

        void fill(T value) {
            ForElements(e)
                *e = value;
        }

        void add(T c) {

            ForElements(e)
                *e += c;
        }

        void add(const BufferTemplate<T>& a) {

            ForLeftRight(a)
                *left += *right;
        }

        void multiply(T scalar) {

            ForElements(e)
                *e *= scalar;
        }

        void multiply(T scalar, BufferTemplate<T>& result) const {
            result.create(numElements);
            auto r = result.begin();
            ForElements(e)
                *r++ = *e * scalar;
        }

        void operator-=(T scalar) {

            ForElements(e)
                *e -= scalar;
        }

        void operator/=(T scalar) {

            ForElements(e)
                *e /= scalar;
        }

        void operator|=(T scalar) {

            ForElements(e)
                *e |= scalar;
        }

        void operator^=(T scalar) {

            ForElements(e)
                *e ^= scalar;
        }

        void operator-=(const BufferTemplate<T>& a) {

            ForLeftRight(a)
                *left -= *right;
        }

        bool operator==(const BufferTemplate<T>& a) const {
            if (numElements != a.numElements)
                return false;
            return memcmp(this->begin(), a.begin(), numElements * sizeof(T)) == 0;
        }

        bool operator!=(const BufferTemplate<T>& a) const {

            if (numElements != a.numElements)
                return true;
            return memcmp(this->begin(), a.begin(), numElements * sizeof(T)) != 0;
        }

        const T*    begin(void) const { return buffer; }
        T*          begin(void) { return buffer; }
        const T*    end(void) const { return endBuffer; }
        T*          end(void) { return endBuffer; }
        size_t      size(void) { return numElements; }

    protected:
        T*          buffer;
        T*          endBuffer;
        size_t      numElements;
};



#pragma mark - ArrayTemplate Definition -

template<typename T>
class ArrayTemplate: public BufferTemplate<T> {

    public:
        ArrayTemplate(void) {

            create(1);  // why?
        }

        
        ArrayTemplate(size_t ne) {

            create(ne);
        }


        ArrayTemplate(const ArrayTemplate<T>& a) {

            BufferTemplate<T>::copy(a);
        }

        ArrayTemplate<T>& operator=(const ArrayTemplate<T>& rhs) {
            if (this != &rhs)
                BufferTemplate<T>::copy(rhs);
            return *this;
        }

        void create(size_t elements) {
            BufferTemplate<T>::create(elements); 
         }

        T operator[](size_t i) const { 
            return this->buffer[i]; 
        }

        T getValue(size_t i) const { return this->buffer[i]; }

        void copy(const ArrayTemplate<T>& other) {
            BufferTemplate<T>::copy(other);
        }
};



#pragma mark - MatrixTemplate Definition -

template<typename T>
class MatrixTemplate : public BufferTemplate<T> {

    public:
        MatrixTemplate(void) : 
            BufferTemplate<T>::BufferTemplate(1) 
        {

            numRows = 1;
            numCols = 1;
        }

        explicit MatrixTemplate(const MatrixTemplate<T>& m) : 
            BufferTemplate<T>::BufferTemplate(m) 
        {

            copy(m);
        }

        explicit MatrixTemplate(size_t nr, size_t nc) : 
            BufferTemplate<T>::BufferTemplate(nr* nc) 
        {

            numRows = nr;
            numCols = nc;
            BufferTemplate<T>::setZero();
        }

        MatrixTemplate<T>& operator=(const MatrixTemplate<T>& rhs) {

            if (this != &rhs)
                copy(rhs);
            return *this;
        }

        T& operator()(size_t r, size_t c) { 
            return this->buffer[r * numCols + c]; 
        }

        const T& operator()(size_t r, size_t c) const { 
            return this->buffer[r * numCols + c]; 
        }

        size_t getNumRows() const {return numRows;}
        size_t getNumCols() const { return numCols; }

        T getValue(size_t r, size_t c) const {

            return this->buffer[r * numCols + c];
        }

        void setValue(size_t r, size_t c, T value) {

            this->buffer[r * numCols + c] = value;
        }

        void print(void) {

            std::cout << "{";
            for (int i = 0; i < numRows; i++)
            {
                std::cout << "{";
                for (int j = 0; j < numCols; j++)
                {
                    std::cout << getValue(i, j);
                    if (j != numCols - 1)
                        std::cout << ",";
                }
                if (i != numRows - 1)
                    std::cout << "}," << std::endl;
            }
            std::cout << "}}" << std::endl;
        }

        void print(std::string title) {

            std::cout << title << std::endl;
            print();
        }

        void copy(const MatrixTemplate<T>& other) {

            create(other.numRows, other.numCols);
            memcpy(this->buffer, other.buffer, this->numElements * sizeof(T));
        }

        void create(size_t nr, size_t nc) {

            BufferTemplate<T>::create(nr * nc);
            numRows = nr;
            numCols = nc;
        }

        void setIdentity() {

            SIZEASSERT(numRows == numCols);
            this->setZero();
            auto cols1 = getNumCols() + 1;
            for (auto end = this->end(), c = this->begin(); c < end; c += cols1)
                *c = 1;
        }

        bool operator==(const MatrixTemplate<T>& m) const {

            if (numRows != m.numRows || numCols != m.numCols)
                return false;
            return BufferTemplate<T>::operator==(m);
        }

        bool operator!=(const MatrixTemplate<T>& m) const {

            if (numRows != m.numRows || numCols != m.numCols)
                return true;
            return BufferTemplate<T>::operator!=(m);
        }

        void transpose(MatrixTemplate<T>& result) {

            result.Create(numCols, numRows);
            auto cell = this->begin();
            for (size_t r = 0; r < numRows; ++r)
                {
                for (int c = 0; c < numCols; ++c)
                    result.setValue(c, r, *cell++);
                }
        }

        void add(const MatrixTemplate<T>& m) {
            SIZEASSERT(numRows == m.numRows && numCols == m.numCols);
            BufferTemplate<T>::add(m);
        }

        void add(const MatrixTemplate<T>& m, MatrixTemplate<T>& result) const {
            SIZEASSERT(numRows == m.numRows && numCols == m.numCols);
            BufferTemplate<T>::add(m, result);
        }

        void multiply(T scalar) {
            BufferTemplate<T>::multiply(scalar);
        }

        void multiply(T scalar, MatrixTemplate<T>& result) const {
            result.create(getNumRows(), getNumCols());
            auto r = result.begin();
            ForElements(e)
                * r++ = *e * scalar;
        }

        void multiply(const MatrixTemplate<T>& m, MatrixTemplate<T>& result) const {

            SIZEASSERT(numCols == m.numRows);
            auto rows = getNumRows();
            auto cols = getNumCols();
            auto mcols = m.getNumCols();
            result.create(rows, mcols);

            auto row = this->begin();
            auto r = result.begin();
            auto mbegin = m.begin();
            for (size_t i = 0; i < rows; i++)
                {
                auto mrow = mbegin;

                for (size_t j = 0; j < mcols; j++)
                    {
                    T sum = 0;
                    auto col = row;
                    auto mcol = mrow;
                    for (size_t k = 0; k < cols; k++)
                        {
                        sum += *col++ * *mcol;
                        mcol += mcols;
                        }

                    *r++ = sum;
                    ++mrow;
                    }

                row += cols;
                }
        }

        void divideByPowerOfTwo(int power) {
            if (power > 0)
              multiply(1.0 / (double)(1 << power));
        }




    protected:
        size_t numRows,
               numCols;
};




#pragma mark - Type Definitions -

typedef ArrayTemplate<int>    IntArray;
typedef ArrayTemplate<double> DoubleArray;

typedef MatrixTemplate<int>    IntMatrix;
typedef MatrixTemplate<double> DoubleMatrix;

#endif
