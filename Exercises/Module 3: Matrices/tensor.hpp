#pragma once

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <vector>

/*
 * This is a modified version of the matrix.h header provided in e-learning
 *
 * It extends matrix calculations to a more easier control over math vectors
 *
 * To create Tensors refer to the TensorBuilder class
 */
template <typename T> class Tensor
{
    template <typename U> friend class TensorBuilder;

  public:
    static Tensor<T> Vector(int n)
    {
        return Tensor<T>(n, 1);
    }

    static Tensor<T> Matrix(int n, int m)
    {
        return Tensor<T>(n, m);
    }

    static Tensor<T> SMatrix(int n)
    {
        return Tensor<T>(n, n);
    }

    static Tensor<T> Identity(int n)
    {
        Tensor<T> out(n, n);
        for (int i = 0; i < n; i++)
            out(i, i) = 1.0;
        return out;
    }

    static Tensor<T> Ones(int n)
    {
        Tensor<T> out(n, n);
        for (int i = 0; i < n; i++)
            out(i, i) = 1.0;
        return out;
    }

    static Tensor<T> One(int n)
    {
        Tensor<T> out(n, 1);
        for (int i = 0; i < n; i++)
            out(i) = 1.0;
        return out;
    }

    static Tensor<T> Random(int n)
    {
        Tensor<T> out(n, n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                out(i, j) = rand() % 100;
        return out;
    }

    static Tensor<T> Random(int n, int m)
    {
        Tensor<T> out(n, m);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                out(i, j) = rand() % 100;
        return out;
    }

    static Tensor<T> FromData(std::vector<T> data, int rows, int cols)
    {
        return Tensor<T>(std::move(data), rows, cols);
    }

    static Tensor<T> FromData(std::vector<T> data)
    {
        return FromData(std::move(data), data.size(), 1);
    }

  private:
    /*
     * Retrieve vector position from tensor coordinates
     * (Starting from 0,0 to nc,nr)
     */
    inline int Site(int i, int j) const
    {
        return i * m_Cols + j;
    }

    Tensor(int rows, int cols)
        : m_Rows(rows), m_Cols(cols), m_Data(rows * cols, 0)
    {
    }

    Tensor(std::vector<T> data, int rows, int cols)
        : m_Rows(rows), m_Cols(cols), m_Data(data)
    {
        assert(data.size() == rows * cols);
    }

  public:
    Tensor() : m_Rows(0), m_Cols(0)
    {
    }

    Tensor(Tensor<T> const &other)
        : m_Rows(other.m_Rows), m_Cols(other.m_Cols), m_Data(other.m_Data)
    {
    }

  public:
    T operator()(int i, int j) const;
    T operator()(int i) const;

    T &operator()(int i, int j);
    T &operator()(int i);

    /* ----------------- COMMON TENSOR OPERATIONS --------------------------- */

    Tensor<T> operator+(Tensor<T> const &b) const;
    Tensor<T> operator-(Tensor<T> const &b) const;

    Tensor<T> Dot(Tensor<T> const &a) const;

    Tensor<T> Dagger();

    T Norm();

    /* ----------------- SOLUTION INVARIANT OPERATIONS ---------------------- */

    void SwapRows(int i, int j);

    /*
     * Linear combination of two rows into a final row
     * of the tensor
     *
     * row_final = row_i + value * row_j
     */
    void LinearCombRows(int i, int j, T value, int final);

    // Action via a scalar
    Tensor<T> operator*(T const &d) const;

    Tensor<T> operator*(T d);

    /* ----------------- SOLUTION INVARIANT OPERATIONS ---------------------- */

    Tensor<T> &operator=(Tensor<T> const &in)
    {
        return *this(m_Data, m_Rows, m_Cols);
    }

    inline bool operator==(Tensor<T> const &b)
    {
        return m_Rows == b.m_Rows && m_Cols == b.m_Cols && m_Data == b.m_Data;
    }

    inline bool &operator==(Tensor<T> const &b) const
    {
        return m_Rows == b.m_Rows && m_Cols == b.m_Cols && m_Data == b.m_Data;
    }

    inline int Cols()
    {
        return m_Cols;
    }

    inline int Rows()
    {
        return m_Rows;
    }

    inline int Cols() const
    {
        return m_Cols;
    }

    inline int Rows() const
    {
        return m_Rows;
    }

    inline std::vector<T> RawData()
    {
        return m_Data;
    }

  private:
    // Values of the Tensor
    std::vector<T> m_Data;

    // Number of rows
    int m_Rows;
    // Number of columns
    int m_Cols;
};

// Implementation of the Tensor class' methods
#include "tensor.inl"
