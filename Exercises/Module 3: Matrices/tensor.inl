#include "tensor.hpp"

template <typename T> void Tensor<T>::SwapRows(int i, int j)
{
    T *tmp = new T[m_Cols];

    // Move slices of the std::vector data on the correct column representation
    // in the std::vector format
    std::copy(tmp, &m_Data[this->Site(i, 0)], sizeof(T) * m_Cols);
    std::copy(&m_Data[this->Site(i, 0)], &m_Data[this->Site(j, 0)],
              sizeof(T) * m_Cols);
    std::copy(&m_Data[this->Site(j, 0)], tmp, sizeof(T) * m_Cols);

    delete[] tmp;
}

template <typename T>
void Tensor<T>::LinearCombRows(int i, int j, T value, int final)
{
    assert(i < m_Rows && j < m_Rows);
    assert(final < m_Rows);

    for (int k = 0; k < m_Cols; k++)
    {
        m_Data[Site(final, k)] =
            m_Data[Site(i, k)] + value * m_Data[Site(j, k)];
    }
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

template <typename T> Tensor<T> Tensor<T>::operator+(Tensor<T> const &b) const
{
    // Tensors should have the same dimension
    assert((m_Rows == b.m_Rows) && (m_Cols == b.m_Cols));

    Tensor<T> out(m_Rows, m_Cols);

    for (int i = 0; i < m_Rows * m_Cols; i++)
        out.m_Data[i] = m_Data[i] + b.m_Data[i];

    return out;
}

template <typename T> Tensor<T> Tensor<T>::Dot(Tensor<T> const &a) const
{
    // Check for compability of the tensors for a dot product
    assert((m_Cols == a.m_Rows));

    Tensor<T> out(m_Rows, a.m_Cols);

    for (int i = 0; i < m_Rows; i++)
        for (int j = 0; j < a.m_Cols; j++)
            for (int k = 0; k < m_Cols; k++)
                out(i, j) += m_Data[this->Site(i, k)] * a(k, j);
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

template <typename T> T Tensor<T>::Norm()
{
    T n = 0.0;

    for (int i = 0; i < m_Rows * m_Cols; i++)
    {
        T h = std::abs(this->m_Data[i]);
        n += h * h;
    }

    return n;
}

template <typename T> Tensor<T> Tensor<T>::operator-(Tensor<T> const &b) const
{
    return *this + (b * (-1.0));
}

template <typename T> Tensor<T> Tensor<T>::operator*(T const &d) const
{
    Tensor<T> out(*this);

    for (int i = 0; i < m_Rows * m_Cols; i++)
        out.m_Data[i] *= d;

    return out;
}

template <typename T> Tensor<T> Tensor<T>::operator*(T d)
{
    Tensor<T> out(*this);
    for (int i = 0; i < m_Rows * m_Cols; i++)
        out.m_Data[i] *= d;

    return out;
}
