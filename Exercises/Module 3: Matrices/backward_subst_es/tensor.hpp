#pragma once

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <vector>

template <typename U> class TensorBuilder;

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

  private:
    /*
     * Retrieve vector position from tensor coordinates
     * (Starting from 0,0 to nc,nr)
     */
    inline int Site(int i, int j) const
    {
        return i * m_Cols + j;
    }

    Tensor(int nrows, int ncols) : m_Rows(nrows), m_Cols(ncols), m_Data(nrows * ncols, 0)
    {
    }

  public:
    Tensor(const Tensor<T> &other) : m_Rows(other.m_Rows), m_Cols(other.m_Cols), m_Data(other.m_Data)
    {
    }

  public:
    T operator()(int i, int j) const;
    T operator()(int i) const;

    T &operator()(int i, int j);
    T &operator()(int i);

    /* Common Tensor operations */

    Tensor<T> operator+(const Tensor<T> &b) const;

    Tensor<T> operator-(const Tensor<T> &b) const
    {
        return *this + (b * (-1.0));
    }

    Tensor<T> Dot(const Tensor<T> &a) const;

    Tensor<T> Dagger();

    double Norm();

    /* Solution invariant operations                   */

    void SwapRows(int i, int j);

    // TODO:
    void SumRows(int i, int j, int final);

    // Action via a scalar
    template <typename K> Tensor<T> operator*(K d) const;

    /****************************************************/

    Tensor<T> &operator=(const Tensor<T> &in)
    {
        m_Rows = in.m_Rows;
        m_Cols = in.m_Cols;
        m_Data = in.m_Data;
        return *this;
    }

    inline bool operator==(const Tensor<T> &b)
    {
        return m_Rows == b.m_Rows && m_Cols == b.m_Cols && m_Data == b.m_Data;
    }

    inline bool &operator==(const Tensor<T> &b) const
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

  private:
    // Values of the Tensor
    std::vector<T> m_Data;

    // Number of rows
    int m_Rows;
    // Number of columns
    int m_Cols;
};

template <typename T> void Tensor<T>::SwapRows(int i, int j)
{
    T *tmp = new T[m_Cols];

    // Move slices of the std::vector data on the correct column representation in
    // the std::vector format
    std::copy(tmp, &m_Data[this->Site(i, 0)], sizeof(T) * m_Cols);
    std::copy(&m_Data[this->Site(i, 0)], &m_Data[this->Site(j, 0)], sizeof(T) * m_Cols);
    std::copy(&m_Data[this->Site(j, 0)], tmp, sizeof(T) * m_Cols);
    delete[] tmp;
}

template <typename T> T Tensor<T>::operator()(int i, int j) const
{
    assert(i < m_Rows && j < m_Cols);

    return m_Data[this->Site(i, j)];
}

template <typename T> T &Tensor<T>::operator()(int i, int j)
{
    assert(i < m_Rows && j < m_Cols);

    return m_Data[this->Site(i, j)];
}

template <typename T> T Tensor<T>::operator()(int i) const
{
    // Tensor should be a vector
    assert(m_Rows == 1 || m_Cols == 1);

    // Check if vector is transposed
    if (m_Cols == 1)
    {
        assert(i < m_Rows);

        return m_Data[this->Site(i, 0)];
    }

    assert(i < m_Cols);

    return m_Data[this->Site(0, i)];
}

template <typename T> T &Tensor<T>::operator()(int i)
{
    // Tensor should be a vector
    assert(m_Rows == 1 || m_Cols == 1);

    // Check if vector is transposed
    if (m_Cols == 1)
    {
        assert(i < m_Rows);

        return m_Data[this->Site(i, 0)];
    }

    assert(i < m_Cols);

    return m_Data[this->Site(0, i)];
}

template <typename T> Tensor<T> Tensor<T>::operator+(const Tensor<T> &b) const
{
    Tensor<T> out(m_Rows, m_Cols);

    // Tensors should have the same dimension
    assert((m_Rows == b.m_Rows) && (m_Cols == b.m_Cols));

    for (int i = 0; i < m_Rows * m_Cols; i++)
        out.m_Data[i] = m_Data[i] + b.m_Data[i];

    return out;
}

template <typename T> Tensor<T> Tensor<T>::Dot(const Tensor<T> &a) const
{
    Tensor<T> out(m_Rows, a.m_Cols);

    // Check for compability of the tensors for a dot product
    assert((m_Cols == a.m_Rows));

    for (int i = 0; i < m_Rows; i++)
        for (int j = 0; j < a.m_Cols; j++)
            for (int k = 0; k < m_Cols; k++)
                out(i, j) += m_Data[this->Site(i, k)] * a(k, j);
    return out;
}

template <typename T> template <typename K> Tensor<T> Tensor<T>::operator*(K d) const
{
    Tensor<T> out(*this);

    for (int i = 0; i < m_Rows * m_Cols; i++)
        out.m_Data[i] *= d;

    return out;
}

template <typename T> Tensor<T> Tensor<T>::Dagger()
{
    Tensor<T> out(m_Cols, m_Rows);

    for (int i = 0; i < m_Cols; i++)
        for (int j = 0; j < m_Rows; j++)
            out(i, j) = conj(m_Data[this->Site(j, i)]);

    return out;
}

template <typename T> double Tensor<T>::Norm()
{
    double n = 0.0;

    for (int i = 0; i < m_Rows * m_Cols; i++)
    {
        double h = std::abs(this->m_Data[i]);
        n += h * h;
    }

    return n;
}
