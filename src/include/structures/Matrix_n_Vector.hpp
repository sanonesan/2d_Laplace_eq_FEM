#pragma once

#ifndef MATRIX_N_VECTOR_HPP
#define MATRIX_N_VECTOR_HPP

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

const double eps = 1e-16;

template <class T>
class Matrix;

template <class T>
class Vector;

//----------------------------------------//
//----------------VECTOR------------------//
//----------------------------------------//

template <class T>
class Vector : public std::vector<T> {
   public:
    friend Matrix<T>;

    Vector<T>();
    explicit Vector<T>(const std::size_t size);
    Vector<T>(const std::vector<T>& vec);
    Vector<T>(const Vector<T>& vec);
    ~Vector<T>(){};

    Vector<T> dot(const Matrix<T>& matrix);

    // vector norms

    T norm_1() const;
    T norm_inf() const;
    T norm_euclid() const;

    // Vector operations

    Vector<T> add(const Vector<T>& another_vector) const;
    Vector<T> subtract(const Vector<T>& another_vector) const;

    T dot(const Vector<T>& another_vector);
    T dot(const std::vector<T>& another_vector);

    // Value operations

    Vector<T> add(const T value) const;
    Vector<T> subtract(const T value) const;
    Vector<T> multiply(const T value) const;
    Vector<T> divide(const T value) const;

    // Vector operations

    Vector<T>& add_this(const Vector<T>& another_vector) const;
    Vector<T>& subtract_this(const Vector<T>& another_vector) const;

    // Value operations

    Vector<T>& add_this(const T value) const;
    Vector<T>& subtract_this(const T value) const;
    Vector<T>& multiply_this(const T value) const;
    Vector<T>& divide_this(const T value) const;

    // Bool operations

    bool operator==(const Vector<T>& another_vector);
    bool operator!=(const Vector<T>& another_vector);

    // Vector operations

    Vector<T> operator+(const Vector<T>& another_vector);
    Vector<T> operator-(const Vector<T>& another_vector);
    T operator*(const Vector<T>& another_vector);

    // Value operations

    Vector<T> operator+(const T value);
    Vector<T> operator-(const T value);
    Vector<T> operator*(const T value);
    Vector<T> operator/(const T value);

    // Matrix operations

    Vector<T>& operator=(const Vector<T>& another_vector);
    Vector<T>& operator+=(const Vector<T>& another_vector);
    Vector<T>& operator-=(const Vector<T>& another_vector);

    // Value operations

    Vector<T>& operator+=(const T value);
    Vector<T>& operator-=(const T value);
    Vector<T>& operator*=(const T value);
    Vector<T>& operator/=(const T value);

    Vector<T>& check_vector_zero();

    friend std::ostream& operator<<(std::ostream& out, const Vector<T>& vec) {
        out << "( " << vec[0];
        for (std::size_t i = 1; i < vec.size() - 1; ++i) {
            out << std::setw(15) << vec[i];
        }
        out << std::setw(15) << vec[vec.size() - 1] << " )^T";

        return out;
    }
};

template <class T>
class Matrix {
   private:
    std::size_t _rows;
    std::size_t _cols;
    std::vector<Vector<T>> _array;

   public:
    friend Vector<T>;

    // Constructors
    Matrix<T>();
    Matrix<T>(std::size_t _n);
    Matrix<T>(std::size_t _rows, std::size_t _cols);
    Matrix<T>(std::vector<Vector<T>> const& _array);
    Matrix<T>(std::vector<std::vector<T>> const& _array);

    // Copy constructor
    Matrix<T>(const Matrix<T>& matrix);
    // Destructor
    ~Matrix<T>(){};

    // getters
    std::size_t get_n() const;
    std::size_t get_rows() const;
    std::size_t get_cols() const;
    std::vector<Vector<T>> get_array() const;

    //-----------RESULT = NEW_MATRIX-----------//

    // Matrix operations

    Matrix<T> add(const Matrix<T>& another_matrix) const;
    Matrix<T> subtract(const Matrix<T>& another_matrix) const;
    Matrix<T> dot(const Matrix<T>& another_matrix) const;
    Matrix<T> transpose() const;
    Matrix<T> inv() const;

    // Value operations

    Matrix<T> add(const T value) const;
    Matrix<T> subtract(const T value) const;
    Matrix<T> multiply(const T value) const;
    Matrix<T> divide(const T value) const;

    //--------INFLUENCE_CURRENT_MATRIX---------//

    // Matrix operations

    Matrix<T>& add_this(const Matrix<T>& another_matrix);
    Matrix<T>& subtract_this(const Matrix<T>& another_matrix);
    Matrix<T>& dot_this(const Matrix<T>& another_matrix);
    Matrix<T>& transpose_this();

    Matrix<T>& make_matrix_identity(const std::size_t n);

    // Value operations

    Matrix<T>& add_this(const T value);
    Matrix<T>& subtract_this(const T value);
    Matrix<T>& multiply_this(const T value);
    Matrix<T>& divide_this(const T value);
    Matrix<T>& nullify();

    // matrix norms

    T norm_1() const;
    T norm_inf() const;

    void read_Square_Matrix(const std::string path);
    void read_Matrix(const std::string path);
    void print();

    // Bool operations

    bool operator==(const Matrix<T>& another_matrix);
    bool operator!=(const Matrix<T>& manother_matrix);

    // Matrix operations

    Matrix<T> operator+(const Matrix<T>& another_matrix);
    Matrix<T> operator-(const Matrix<T>& another_matrix);
    Matrix<T> operator*(const Matrix<T>& another_matrix);

    // Value operations

    Matrix<T> operator+(const T value);
    Matrix<T> operator-(const T value);
    Matrix<T> operator*(const T value);

    // Matrix operations

    Matrix<T>& operator=(const Matrix<T>& another_matrix);
    Matrix<T>& operator+=(const Matrix<T>& another_matrix);
    Matrix<T>& operator-=(const Matrix<T>& another_matrix);
    Matrix<T>& operator*=(const Matrix<T>& another_matrix);

    // Value operations

    Matrix<T>& operator+=(const T value);
    Matrix<T>& operator-=(const T value);
    Matrix<T>& operator*=(const T value);

    Matrix<T>& check_matrix_zero();

    Vector<T>& operator[](const std::size_t i);

    Vector<T> dot(const Vector<T>& vector);
    Vector<T> dot(const std::vector<T>& vector);

    template <typename T1>
    friend void read_System(const std::string path, Matrix<T1>& matrix, Vector<T1>& vector);

    friend std::ostream& operator<<(std::ostream& out, Matrix<T>& matrix) {
        if (matrix.get_rows() == 0 || matrix.get_cols() == 0)
            throw std::invalid_argument("null matrix");

        for (std::size_t i = 0; i < matrix.get_rows(); ++i) {
            for (std::size_t j = 0; j < matrix.get_cols(); ++j)
                out << std::setw(15) << matrix.get_array()[i][j];
            out << "\n";
        }
        out << "\n";

        return out;
    };
};

#endif

//----------------------------------------//
//----------------MATRIX------------------//
//----------------------------------------//

template <class T>
Matrix<T>::Matrix() {
    this->_rows = 0;
    this->_cols = 0;
}

template <class T>
Matrix<T>::Matrix(std::size_t _n) {
    this->_cols = _n;
    this->_rows = _n;
    this->_array = std::vector<Vector<T>>(_n, std::vector<T>(_n));
}

template <class T>
Matrix<T>::Matrix(std::size_t _rows, std::size_t _cols) {
    this->_rows = _rows;
    this->_cols = _cols;
    this->_array = std::vector<Vector<T>>(_rows, std::vector<T>(_cols));

    for (std::size_t i = 0; i < this->_rows; ++i) {
        for (std::size_t j = 0; j < this->_cols; ++j) {
            this->_array[i][j] = 0.;
        }
    }
}

template <class T>
Matrix<T>::Matrix(std::vector<Vector<T>> const& _array) {
    if (_array.size() == 0 || _array[0].size() == 0)
        throw std::invalid_argument("Size of _array should be > 0");

    this->_rows = _array.size();
    this->_cols = _array[0].size();
    this->_array = _array;
}

template <class T>
Matrix<T>::Matrix(std::vector<std::vector<T>> const& _array) {
    if (_array.size() == 0 || _array[0].size() == 0)
        throw std::invalid_argument("Size of _array should be > 0");

    std::vector<Vector<T>> res;
    this->_rows = _array.size();
    this->_cols = _array[0].size();

    this->_array.resize(this->_rows);
    for (std::size_t i = 0; i < this->_rows; ++i)
        this->_array[i] = _array[i];
}

// Copy constructor
template <class T>
Matrix<T>::Matrix(const Matrix<T>& matrix) {
    this->_rows = matrix._rows;
    this->_cols = matrix._cols;
    this->_array = matrix._array;
}

// Getters

template <class T>
std::size_t Matrix<T>::get_rows() const {
    if (this->_rows == 0)
        std::cout << "null matrix!";

    return this->_rows;
}

template <class T>
std::size_t Matrix<T>::get_cols() const {
    if (this->_cols == 0)
        std::cout << "null matrix!";
    return this->_cols;
}

template <class T>
std::vector<Vector<T>> Matrix<T>::get_array() const {
    if (this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    return this->_array;
}

/*
/////////////////////////////////////////////
//-----------RESULT = NEW_MATRIX-----------//
/////////////////////////////////////////////
*/

// Matrix operations

template <class T>
Matrix<T> Matrix<T>::add(const Matrix<T>& another_matrix) const {
    if (this->_rows == 0 || this->_cols == 0 || another_matrix._rows == 0 || another_matrix._cols == 0)
        throw std::invalid_argument("null matrix");

    if (this->_rows != another_matrix._rows && this->_cols != another_matrix._cols)
        throw std::invalid_argument("These matrixs cannot be added first._cols != second._rows && first._rows != second._cols");

    Matrix<T> result_matrix(this->_array);

    T tmp;
    for (std::size_t i = 0; i < this->_rows; ++i) {
        for (std::size_t j = 0; j < this->_cols; ++j) {
            tmp = result_matrix._array[i][j];
            tmp += another_matrix._array[i][j];

            if (fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;

            result_matrix._array[i][j] = tmp;
        }
    }

    return result_matrix;
}

template <class T>
Matrix<T> Matrix<T>::subtract(const Matrix<T>& another_matrix) const {
    if (this->_rows == 0 || this->_cols == 0 || another_matrix._rows == 0 || another_matrix._cols == 0)
        throw std::invalid_argument("null matrix");

    if (this->_rows != another_matrix._rows && this->_cols != another_matrix._cols)
        throw std::invalid_argument("These matrixs cannot be substructed first._cols != second._rows && first._rows != second._cols");

    Matrix<T> result_matrix(this->_array);

    T tmp;
    for (std::size_t i = 0; i < this->_rows; ++i) {
        for (std::size_t j = 0; j < this->_cols; ++j) {
            tmp = result_matrix._array[i][j];
            tmp -= another_matrix._array[i][j];

            if (fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;

            result_matrix._array[i][j] = tmp;
        }
    }

    return result_matrix;
}

template <class T>
Matrix<T> Matrix<T>::dot(const Matrix<T>& another_matrix) const {
    if (this->_rows == 0 || this->_cols == 0 || another_matrix._rows == 0 || another_matrix._cols == 0)
        throw std::invalid_argument("null matrix");

    if (this->_cols != another_matrix._rows) {
        throw std::invalid_argument("These matrixs cannot be multiplied first._cols != second._rows");
    }

    if (this->_cols == 0 || this->_rows == 0 || another_matrix._rows == 0 || another_matrix._cols == 0) {
        throw std::invalid_argument("These matrixs cannot be multiplied null matrix");
    }

    Matrix<T> result_matrix(this->_rows, another_matrix._cols);

    T tmp = 0;
    for (std::size_t i = 0; i < result_matrix._rows; ++i) {
        for (std::size_t j = 0; j < result_matrix._cols; ++j) {
            tmp = 0;

            for (std::size_t k = 0; k < result_matrix._cols; ++k)
                tmp += this->_array[i][k] * another_matrix._array[k][j];

            if (fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;

            result_matrix._array[i][j] = tmp;
        }
    }

    return result_matrix;
}

template <class T>
Matrix<T> Matrix<T>::transpose() const {
    if (this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    Matrix<T> result_matrix(this->_cols, this->_rows);

    for (std::size_t i = 0; i < this->_rows; ++i)
        for (std::size_t j = i; j < this->_cols; ++j) {
            result_matrix._array[i][j] = this->_array[j][i];
            result_matrix._array[j][i] = this->_array[i][j];
        }

    return result_matrix;
}

// Value operations
template <class T>
Matrix<T> Matrix<T>::add(const T value) const {
    if (this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    Matrix<T> result_matrix(this->_rows, this->_cols);

    T tmp;
    for (std::size_t i = 0; i < this->_rows; ++i) {
        for (std::size_t j = 0; j < this->_cols; ++j) {
            tmp = this->_array[i][j] + value;

            if (fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;

            result_matrix._array[i][j] = tmp;
        }
    }

    return result_matrix;
}

template <class T>
Matrix<T> Matrix<T>::subtract(const T value) const {
    if (this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    Matrix<T> result_matrix(this->_rows, this->_cols);

    T tmp;
    for (std::size_t i = 0; i < this->_rows; ++i) {
        for (std::size_t j = 0; j < this->_cols; ++j) {
            tmp = this->_array[i][j] - value;

            if (fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;

            result_matrix._array[i][j] = tmp;
        }
    }

    return result_matrix;
}

template <class T>
Matrix<T> Matrix<T>::multiply(const T value) const {
    if (this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    Matrix<T> result_matrix(this->_rows, this->_cols);

    T tmp;
    for (std::size_t i = 0; i < this->_rows; ++i) {
        for (std::size_t j = 0; j < this->_cols; ++j) {
            tmp = this->_array[i][j] * value;

            if (fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;

            result_matrix._array[i][j] = tmp;
        }
    }

    return result_matrix;
}

template <class T>
Matrix<T> Matrix<T>::divide(const T value) const {
    if (this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    Matrix<T> result_matrix(this->_rows, this->_cols);

    T tmp;
    for (std::size_t i = 0; i < this->_rows; ++i) {
        for (std::size_t j = 0; j < this->_cols; ++j) {
            tmp = this->_array[i][j] / value;

            if (fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;

            result_matrix._array[i][j] = tmp;
        }
    }

    return result_matrix;
}

/*
/////////////////////////////////////////////
//--------INFLUENCE_CURRENT_MATRIX---------//
/////////////////////////////////////////////
*/

// Matrix operations

template <class T>
Matrix<T>& Matrix<T>::add_this(const Matrix<T>& another_matrix) {
    if (this->_rows == 0 || this->_cols == 0 || another_matrix._rows == 0 || another_matrix._cols == 0)
        throw std::invalid_argument("null matrix");

    if (this->_rows != another_matrix._rows && this->_cols != another_matrix._cols)
        throw std::invalid_argument("These matrixs cannot be added first._cols != second._rows && first._rows != second._cols");

    T tmp;
    for (std::size_t i = 0; i < this->_rows; ++i) {
        for (std::size_t j = 0; j < this->_cols; ++j) {
            tmp = this->_array[i][j];
            tmp += another_matrix._array[i][j];

            if (fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;

            this->_array[i][j] = tmp;
        }
    }

    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::subtract_this(const Matrix<T>& another_matrix) {
    if (this->_rows == 0 || this->_cols == 0 || another_matrix._rows == 0 || another_matrix._cols == 0)
        throw std::invalid_argument("null matrix");

    if (this->_rows != another_matrix._rows && this->_cols != another_matrix._cols)
        throw std::invalid_argument("These matrixs cannot be substructed first._cols != second._rows && first._rows != second._cols");

    T tmp;
    for (std::size_t i = 0; i < this->_rows; ++i) {
        for (std::size_t j = 0; j < this->_cols; ++j) {
            tmp = this->_array[i][j];
            tmp -= another_matrix._array[i][j];

            if (fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;

            this->_array[i][j] = tmp;
        }
    }

    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::dot_this(const Matrix<T>& another_matrix) {
    if (this->_rows == 0 || this->_cols == 0 || another_matrix._rows == 0 || another_matrix._cols == 0)
        throw std::invalid_argument("null matrix");

    if (this->_cols != another_matrix._rows) {
        throw std::invalid_argument("These matrixs cannot be multiplied first._cols != second._rows");
    }

    if (this->_cols == 0 || this->_rows == 0 || another_matrix._rows == 0 || another_matrix._cols == 0) {
        throw std::invalid_argument("These matrixs cannot be multiplied null matrix");
    }

    std::vector<Vector<T>> result_array(this->_rows, std::vector<T>(another_matrix._cols));

    T tmp = 0;
    for (std::size_t i = 0; i < this->_rows; ++i) {
        for (std::size_t j = 0; j < another_matrix._cols; ++j) {
            tmp = 0.;

            for (std::size_t k = 0; k < another_matrix._cols; ++k)
                tmp += this->_array[i][k] * another_matrix._array[k][j];

            if (fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;

            result_array[i][j] = tmp;
        }
    }

    this->_cols = another_matrix._cols;
    this->_array = result_array;

    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::transpose_this() {
    if (this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    if (this->_rows != this->_cols) {
        Matrix result_matrix(this->_cols, this->_rows);

        std::vector<Vector<T>> result_array(this->_cols, std::vector<T>(this->_rows));

        for (std::size_t i = 0; i < this->_rows; ++i)
            for (std::size_t j = i; j < this->_cols; ++j) {
                result_array[j][i] = this->_array[i][j];
                result_array[i][j] = this->_array[j][i];
            }
        this->_rows = result_array.size();
        this->_cols = result_array[0].size();
        this->_array = result_array;

    } else {
        T tmp;
        for (std::size_t i = 0; i < this->_rows; ++i)
            for (std::size_t j = i + 1; j < this->_cols; ++j) {
                tmp = this->_array[i][j];
                this->_array[i][j] = this->_array[j][i];
                this->_array[j][i] = tmp;
            }
    }

    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::make_matrix_identity(const std::size_t n) {
    if (n == 0)
        throw std::invalid_argument("n must be > 0");

    this->_rows = n;
    this->_cols = n;

    this->_array.resize(n);

    for (std::size_t i = 0; i < n; ++i) {
        this->_array[i].resize(n);
        for (std::size_t j = 0; j < n; ++j) {
            this->_array[i][j] = (i == j) ? 1. : 0.;
        }
    }

    return (*this);
}

// Value operations

template <class T>
Matrix<T>& Matrix<T>::add_this(const T value) {
    if (this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    T tmp;

    for (std::size_t i = 0; i < this->_rows; ++i) {
        for (std::size_t j = 0; j < this->_cols; ++j) {
            tmp = this->_array[i][j] + value;
            if (fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;

            this->_array[i][j] = tmp;
        }
    }

    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::subtract_this(const T value) {
    if (this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    T tmp;
    for (std::size_t i = 0; i < this->_rows; ++i) {
        for (std::size_t j = 0; j < this->_cols; ++j) {
            tmp = this->_array[i][j] - value;

            if (fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;

            this->_array[i][j] = tmp;
        }
    }

    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::multiply_this(const T value) {
    if (this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    T tmp;
    for (std::size_t i = 0; i < this->_rows; ++i) {
        for (std::size_t j = 0; j < this->_cols; ++j) {
            tmp = this->_array[i][j] * value;

            if (fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;

            this->_array[i][j] = tmp;
        }
    }

    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::divide_this(const T value) {
    if (this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    T tmp;
    for (std::size_t i = 0; i < this->_rows; ++i) {
        for (std::size_t j = 0; j < this->_cols; ++j) {
            tmp = this->_array[i][j] / value;

            if (fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;

            this->_array[i][j] = tmp;
        }
    }

    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::nullify() {

    if (this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    std::vector<T> tmp;

    for(std::size_t i = 0; i < this->_cols; ++i)
        tmp.push_back(0.);

    for (std::size_t i = 0; i < this->_rows; ++i)
        this->_array[i].assign(tmp.begin(), tmp.end());
    

    return *this;
}

template <class T>
T Matrix<T>::norm_1() const {
    T norm = 0.;
    T tmp = 0.;
    for (std::size_t i = 0; i < this->_rows; ++i) {
        tmp = 0.;
        for (std::size_t j = 0; j < this->_cols; ++j) {
            tmp += fabs(this->_array[j][i]);
        }
        if (norm < tmp) norm = tmp;
    }
    if (fabs(norm) < eps && norm != 0.) {
        norm = 0.;
        if (1 / norm < 0) norm *= -1;
    }
    return norm;
}

template <class T>
T Matrix<T>::norm_inf() const {
    T norm = 0.;
    T tmp = 0.;
    for (std::size_t i = 0; i < this->_rows; ++i) {
        tmp = 0.;
        for (std::size_t j = 0; j < this->_cols; ++j) {
            tmp += fabs(this->_array[i][j]);
        }
        if (norm < tmp) norm = tmp;
    }
    if (fabs(norm) < eps && norm != 0.) {
        norm = 0.;
        if (1 / norm < 0) norm *= -1;
    }
    return norm;
}

template <class T>
void Matrix<T>::read_Square_Matrix(const std::string path) {
    std::ifstream fin;
    fin.open(path);
    if (!fin.is_open()) {
        throw std::invalid_argument("Smth's wrong with path");
    }
    fin >> this->_rows;
    this->_cols = this->_rows;
    fin.get();
    std::vector<Vector<T>> _array(this->_rows, std::vector<T>(this->_cols));

    for (std::size_t i = 0; i < this->_rows; ++i) {
        for (std::size_t j = 0; j < this->_cols; ++j) {
            fin >> _array[i][j];
        }
    }

    fin.close();
    this->_array = _array;
}

template <class T>
void Matrix<T>::read_Matrix(const std::string path) {
    std::ifstream fin;
    fin.open(path);
    if (!fin.is_open()) {
        throw std::invalid_argument("Smth's wrong with path");
    }
    fin >> this->_rows;
    fin >> this->_cols;
    fin.get();
    std::vector<Vector<T>> colsarray(this->_rows, std::vector<T>(this->_cols));

    for (std::size_t i = 0; i < this->_rows; ++i) {
        for (std::size_t j = 0; j < this->_cols; ++j) {
            fin >> _array[i][j];
        }
    }

    fin.close();
    this->_array = _array;
}

template <typename T>
void read_System(const std::string path, Matrix<T>& matrix, Vector<T>& vector) {
    std::ifstream fin;
    fin.open(path);
    if (!fin.is_open()) {
        throw std::invalid_argument("Smth's wrong with path");
    }

    fin >> matrix._rows;
    matrix._cols = matrix._rows;

    fin.get();

    matrix._array.resize(matrix._rows);
    matrix._cols = matrix._rows;

    vector.resize(matrix._rows);

    std::vector<Vector<T>> _array(matrix._rows, std::vector<T>(matrix._cols));

    for (std::size_t i = 0; i < matrix._rows; ++i) {
        for (std::size_t j = 0; j < matrix._cols; ++j) {
            fin >> _array[i][j];
        }
        fin >> vector[i];
    }

    fin.close();
    matrix._array = _array;

    return;
}

template <class T>
Matrix<T>& Matrix<T>::check_matrix_zero() {
    for (std::size_t i = 0; i < (*this)._array.size(); ++i) {
        (*this)._array[i].check_vector_zero();
    }

    return (*this);
}

template <class T>
void Matrix<T>::print() {
    if (this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    for (std::size_t i = 0; i < this->_rows; ++i) {
        for (std::size_t j = 0; j < this->_cols; ++j)
            std::cout << std::setw(15) << this->_array[i][j];
        std::cout << "\n";
    }
    std::cout << "\n";
}

template <class T>
bool Matrix<T>::operator==(const Matrix<T>& another_matrix) {
    if (this->_rows == 0 || this->_cols == 0 || another_matrix._rows == 0 || another_matrix._cols == 0)
        return false;

    if (this->_rows != another_matrix._rows || this->_cols != another_matrix._cols)
        return false;

    for (std::size_t i = 0; i < this->_rows; ++i) {
        for (std::size_t j = 0; j < this->_cols; ++j) {
            if (fabs(this->_array[i][j] - another_matrix._array[i][j]) > eps)
                return false;
        }
    }

    return true;
}

template <class T>
bool Matrix<T>::operator!=(const Matrix<T>& another_matrix) {
    return !operator==(another_matrix);
}

// Matrix operations;

template <class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& another_matrix) {
    return this->add(another_matrix);
}

template <class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& another_matrix) {
    return this->subtract(another_matrix);
}

template <class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& another_matrix) {
    return this->dot(another_matrix);
}

// Value operations

template <class T>
Matrix<T> Matrix<T>::operator+(const T value) {
    return this->add(value);
}

template <class T>
Matrix<T> Matrix<T>::operator-(const T value) {
    return this->subtract(value);
}

template <class T>
Matrix<T> Matrix<T>::operator*(const T value) {
    return this->multiply(value);
}

// Matrix operations;
template <class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& another_matrix) {
    this->_rows = another_matrix._rows;
    this->_cols = another_matrix._cols;
    this->_array = another_matrix._array;
    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& another_matrix) {
    return this->add_this(another_matrix);
}

template <class T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& another_matrix) {
    return this->subtract_this(another_matrix);
}

template <class T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& another_matrix) {
    return this->dot_this(another_matrix);
}

// Value operations

template <class T>
Matrix<T>& Matrix<T>::operator+=(const T value) {
    return this->add_this(value);
}

template <class T>
Matrix<T>& Matrix<T>::operator-=(const T value) {
    return this->subtract_this(value);
}

template <class T>
Matrix<T>& Matrix<T>::operator*=(const T value) {
    return this->multiply_this(value);
}

template <class T>
Vector<T>& Matrix<T>::operator[](const std::size_t i) {
    if (i >= this->_rows)
        throw std::invalid_argument("index out of range");

    return this->_array[i];
};

template <class T>
Vector<T> Matrix<T>::dot(const Vector<T>& vector) {
    if (this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    if (vector.size() == 0)
        throw std::invalid_argument("null vector");

    if (this->_cols != vector.size())
        throw std::invalid_argument("matrix.cols != vector.size");

    Vector<T> result_vector(this->_rows);

    T tmp;

    for (std::size_t i = 0; i < this->_rows; ++i) {
        tmp = 0.;
        for (std::size_t j = 0; j < this->_cols; ++j) {
            tmp += this->_array[i][j] * vector[j];
        }
        if (fabs(tmp) < eps)
            if (1 / tmp < 0) tmp *= -1;

        result_vector[i] = tmp;
    }

    return result_vector;
}

template <class T>
Vector<T> Matrix<T>::dot(const std::vector<T>& vector) {
    Vector<T> tmp(vector);
    return this->dot(tmp);
};

//----------------------------------------//
//----------------VECTOR------------------//
//----------------------------------------//

template <class T>
Vector<T>::Vector() : std::vector<T>(){};

template <class T>
Vector<T>::Vector(const std::size_t size) : std::vector<T>() {
    (*this).resize(size);

    for (std::size_t i = 0; i < size; ++i)
        (*this)[i] = 0.;
};

template <class T>
Vector<T>::Vector(const std::vector<T>& vec) {
    (*this).resize(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i)
        (*this)[i] = vec[i];
};

template <class T>
Vector<T>::Vector(const Vector<T>& vec) {
    (*this).resize(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i)
        (*this)[i] = vec[i];
};

template <class T>
Vector<T> Vector<T>::dot(const Matrix<T>& matrix) {
    if (matrix._rows == 0 || matrix._cols == 0)
        throw std::invalid_argument("null matrix");

    if ((*this).size() == 0)
        throw std::invalid_argument("null vector");

    if ((*this).size() != matrix._rows)
        throw std::invalid_argument("vector._size != matrix._rows");

    Vector<T> result_vector((*this).size());

    T tmp;

    for (std::size_t i = 0; i < (*this).size(); ++i) {
        tmp = 0.;
        for (std::size_t j = 0; j < matrix._rows; ++j) {
            tmp += (*this)[i] * matrix._array[i][j];
        }
        if (fabs(tmp) < eps)
            if (1 / tmp < 0) tmp *= -1;

        result_vector[i] = tmp;
    }

    return result_vector;
}

/*
        Октаэдрическая норма ( || * ||_{1})
*/
template <class T>
T Vector<T>::norm_1() const {
    T norm = 0;
    for (std::size_t i = 0; i < (*this).size(); ++i)
        norm += fabs((*this)[i]);

    return norm < eps ? 0. : norm;
}

/*
        Кубическая норма ( || * ||_{inf})
*/
template <class T>
T Vector<T>::norm_inf() const {
    T norm = fabs((*this)[0]);
    for (std::size_t i = 1; i < (*this).size(); ++i)
        if (fabs(norm - fabs((*this)[i]) < eps))
            norm = fabs((*this)[i]);

    return norm < eps ? 0. : norm;
}

/*
        Евклидова норма ( || * ||_{eucl})
*/
template <class T>
T Vector<T>::norm_euclid() const {
    T norm = 0;
    for (std::size_t i = 0; i < (*this).size(); ++i)
        norm += fabs((*this)[i] * (*this)[i]);

    return norm < eps ? 0. : sqrt(norm);
}

// Bool operations
template <class T>
bool Vector<T>::operator==(const Vector<T>& another_vector) {
    if ((*this).size() == 0 && another_vector.size() == 0)
        return true;

    if ((*this).size() != another_vector.size())
        return false;

    for (std::size_t i = 0; i < (*this).size(); ++i) {
        if ((*this)[i] != another_vector[i])
            return false;
    }

    return true;
};

template <class T>
bool Vector<T>::operator!=(const Vector<T>& another_vector) {
    return !((*this) == another_vector);
};

// Matrix operations

template <class T>
Vector<T> Vector<T>::operator+(const Vector<T>& another_vector) {
    if ((*this).size() == 0 && another_vector.size() == 0)
        return (*this);

    if ((*this).size() != another_vector.size())
        throw std::invalid_argument(" different size");

    Vector<T> result_vector((*this).size());

    for (std::size_t i = 0; i < (*this).size(); ++i) {
        result_vector[i] = (*this)[i] + another_vector[i];
    }

    return result_vector;
};

template <class T>
Vector<T> Vector<T>::operator-(const Vector<T>& another_vector) {
    if ((*this).size() == 0 && another_vector.size() == 0)
        return (*this);

    if ((*this).size() != another_vector.size())
        throw std::invalid_argument(" different size");

    Vector<T> result_vector((*this).size());

    for (std::size_t i = 0; i < (*this).size(); ++i) {
        result_vector[i] = (*this)[i] - another_vector[i];
    }

    return result_vector;
};

template <class T>
T Vector<T>::operator*(const Vector<T>& another_vector) {
    if ((*this).size() == 0 || another_vector.size() == 0)
        throw std::invalid_argument("null vector");

    if ((*this).size() != another_vector.size())
        throw std::invalid_argument(" vectors have different size");

    T tmp = 0;
    for (std::size_t i = 0; i < (*this).size(); ++i) {
        tmp += (*this)[i] * another_vector[i];
    }

    return tmp;
};

// //Value operations
template <class T>
Vector<T> Vector<T>::operator+(const T value) {
    if ((*this).size() == 0)
        throw std::invalid_argument(" different size");

    Vector<T> result_vector((*this).size());

    for (std::size_t i = 0; i < (*this).size(); ++i) {
        result_vector[i] = (*this)[i] + value;
    }

    return result_vector;
};

template <class T>
Vector<T> Vector<T>::operator-(const T value) {
    if ((*this).size() == 0)
        throw std::invalid_argument(" different size");

    Vector<T> result_vector((*this).size());

    for (std::size_t i = 0; i < (*this).size(); ++i) {
        result_vector[i] = (*this)[i] - value;
    }

    return result_vector;
};

template <class T>
Vector<T> Vector<T>::operator*(const T value) {
    if ((*this).size() == 0)
        throw std::invalid_argument("null vector");

    Vector<T> result_vector((*this).size());

    for (std::size_t i = 0; i < (*this).size(); ++i) {
        result_vector[i] = (*this)[i] * value;
    }

    return result_vector;
};

template <class T>
Vector<T> Vector<T>::operator/(const T value) {
    if ((*this).size() == 0)
        throw std::invalid_argument("null vector");

    Vector<T> result_vector((*this).size());

    for (std::size_t i = 0; i < (*this).size(); ++i) {
        result_vector[i] = (*this)[i] / value;
    }

    return result_vector;
};

// Vector operations

template <class T>
Vector<T>& Vector<T>::operator=(const Vector<T>& another_vector) {
    if ((*this).size() == 0 && another_vector.size() == 0)
        return (*this);

    (*this).resize(another_vector.size());
    (*this).assign(another_vector.begin(), another_vector.end());

    return (*this);
};

template <class T>
Vector<T>& Vector<T>::operator+=(const Vector<T>& another_vector) {
    if ((*this).size() == 0 && another_vector.size() == 0)
        return (*this);

    if ((*this).size() != another_vector.size())
        throw std::invalid_argument(" different size");

    for (std::size_t i = 0; i < (*this).size(); ++i) {
        (*this)[i] += another_vector[i];
    }

    return (*this);
};

template <class T>
Vector<T>& Vector<T>::operator-=(const Vector<T>& another_vector) {
    if ((*this).size() == 0 && another_vector.size() == 0)
        return (*this);

    if ((*this).size() != another_vector.size())
        throw std::invalid_argument(" different size");

    for (std::size_t i = 0; i < (*this).size(); ++i) {
        (*this)[i] -= another_vector[i];
    }

    return (*this);
};

// Value operations

template <class T>
Vector<T>& Vector<T>::operator+=(const T value) {
    if ((*this).size() == 0)
        throw std::invalid_argument("null vector");

    for (std::size_t i = 0; i < (*this).size(); ++i) {
        (*this)[i] += value;
    }

    return (*this);
};

template <class T>
Vector<T>& Vector<T>::operator-=(const T value) {
    if ((*this).size() == 0)
        throw std::invalid_argument("null vector");

    for (std::size_t i = 0; i < (*this).size(); ++i) {
        (*this)[i] -= value;
    }

    return (*this);
};

template <class T>
Vector<T>& Vector<T>::operator*=(const T value) {
    if ((*this).size() == 0)
        throw std::invalid_argument("null vector");

    for (std::size_t i = 0; i < (*this).size(); ++i) {
        (*this)[i] *= value;
    }

    return (*this);
};

template <class T>
Vector<T>& Vector<T>::operator/=(const T value) {
    if ((*this).size() == 0)
        throw std::invalid_argument("null vector");

    for (std::size_t i = 0; i < (*this).size(); ++i) {
        (*this)[i] /= value;
    }

    return (*this);
};

template <class T>
T Vector<T>::dot(const Vector<T>& another_vector) {
    return (*this) * another_vector;
};

template <class T>
T Vector<T>::dot(const std::vector<T>& another_vector) {
    Vector<T> tmp(another_vector);

    return (*this) * tmp;
};

// Vector operations

template <class T>
Vector<T> Vector<T>::add(const Vector<T>& another_vector) const {
    return (*this) + another_vector;
};

template <class T>
Vector<T> Vector<T>::subtract(const Vector<T>& another_vector) const {
    return (*this) - another_vector;
};

// Value operations

template <class T>
Vector<T> Vector<T>::add(const T value) const {
    return (*this) + value;
};

template <class T>
Vector<T> Vector<T>::subtract(const T value) const {
    return (*this) - value;
};

template <class T>
Vector<T> Vector<T>::multiply(const T value) const {
    return (*this) * value;
};

template <class T>
Vector<T> Vector<T>::divide(const T value) const {
    return (*this) / value;
};

// Vector operations

template <class T>
Vector<T>& Vector<T>::add_this(const Vector<T>& another_vector) const {
    (*this) += another_vector;
    return (*this);
};

template <class T>
Vector<T>& Vector<T>::subtract_this(const Vector<T>& another_vector) const {
    (*this) -= another_vector;
    return (*this);
};

// Value operations

template <class T>
Vector<T>& Vector<T>::add_this(const T value) const {
    (*this) += value;
    return (*this);
};

template <class T>
Vector<T>& Vector<T>::subtract_this(const T value) const {
    (*this) -= value;
    return (*this);
};

template <class T>
Vector<T>& Vector<T>::multiply_this(const T value) const {
    (*this) *= value;
    return (*this);
};

template <class T>
Vector<T>& Vector<T>::divide_this(const T value) const {
    (*this) /= value;
    return (*this);
};

template <class T>
Vector<T>& Vector<T>::check_vector_zero() {
    const T _eps = 1e-14;
    for (std::size_t i = 0; i < (*this).size(); ++i) {
        if (fabs((*this)[i]) < _eps) {
            (*this)[i] = 0.;
            (*this)[i] *= (1 / (*this)[i] > 0) ? (1) : (-1);
        }
    }

    return (*this);
};

template <class T>
Matrix<T> Matrix<T>::inv() const {
    auto reverse_course = [](
                        Matrix<T>& A, 
                        Vector<T>& b, 
                        Vector<T>& solution) -> int {
        solution = b;
        std::size_t n = solution.size();

        solution[n - 1] /= A[n - 1][n - 1];

        for (std::size_t i = 0; i <= n - 2; ++i) {
            for (std::size_t j = n - 2 - i + 1; j < n; ++j) {
                solution[n - 2 - i] -= solution[j] * A[n - 2 - i][j];
            }
            solution[n - 2 - i] /= A[n - 2 - i][n - 2 - i];
        }

        solution.check_vector_zero();

        return 0;
    };

    auto coefs = [](
                    const std::size_t& k,
                    const std::size_t& l,
                    T& c,
                    T& s,
                    Matrix<T>& A) -> int {

        std::size_t n = A.get_rows();

        if (k < n && l < n) {
            T temp = sqrt(pow(A[k][k], 2) + pow(A[l][k], 2));

            c = A[k][k] / temp;

            s = A[l][k] / temp;

            if (fabs(c) < eps)
                c = (1 / c > 0) ? c : c * (-1);
            if (fabs(s) < eps)
                s = (1 / s > 0) ? s : s * (-1);
        }

        return 0;
    };

    auto QR_decomposion_method = [coefs, reverse_course](

                                    const Matrix<T>& A,
                                    Matrix<T>& Q,
                                    Matrix<T>& R) -> int {
        std::size_t n = A.get_rows();

        T c = 0., s = 0.;

        Q.make_matrix_identity(n);
        R = A;

        T a = 0., b = 0.;

        for (std::size_t i = 0; i < n - 1; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                coefs(i, j, c, s, R);

                for (std::size_t k = 0; k < n; ++k) {
                    a = R[i][k];
                    b = R[j][k];
                    R[i][k] = (c * a + s * b);
                    R[j][k] = (-s * a + c * b);

                    a = Q[i][k];
                    b = Q[j][k];
                    Q[i][k] = (c * a + s * b);
                    Q[j][k] = (-s * a + c * b);
                }
            }
        }

        Q.check_matrix_zero();
        R.check_matrix_zero();

        return 0;
    };

    Matrix<T> Q, R;

    QR_decomposion_method((*this), Q, R);

    Matrix<T> InvMatrix((*this));

    std::size_t n = this->get_rows();
    Vector<T> b(n);
    Vector<T> tmp(n);
    Vector<T> tmp_solution(n);
    // Matrix<T> B(Q);

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            b[j] = (i == j) ? 1. : 0.;
        }

        tmp = Q.dot(b);

        reverse_course(R, tmp, InvMatrix[i]);
    }
    InvMatrix.transpose_this();
    InvMatrix.check_matrix_zero();

    return InvMatrix;
}
