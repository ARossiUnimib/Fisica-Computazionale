#pragma once

#include <cassert>
#include <cmath>
#include <cwchar>
#include <vector>

#include "../utils.hpp"

/**
 * @brief Class to store data of a function
 *
 * @tparam T precision
 */
template <typename T> class FunctionData;

namespace FunctionUtils
{

/**
 * @brief Create a polynomial function
 *
 * @tparam T
 * @param coefficients
 * @param x values
 * @return FunctionData<T>
 */
template <typename T> FunctionData<T> Polynomial(std::vector<T> coefficients, const Utils::Range<T> &x);

} // namespace FunctionUtils

template <typename T> class FunctionData
{
  public:
    FunctionData() : m_xValues(), m_fValues()
    {
    }

    FunctionData(std::size_t size)
    {
        m_xValues.reserve(size);
        m_fValues.reserve(size);
    }

    FunctionData(std::vector<T> xValues, std::vector<T> fValues)
        : m_xValues(std::move(xValues)), m_fValues(std::move(fValues))
    {
        assert(m_xValues.size() == m_fValues.size());
    }

    std::size_t Size() const
    {
        return m_xValues.size();
    }

    T X(std::size_t index) const
    {
        return m_xValues[index];
    }

    T F(std::size_t index) const
    {
        return m_fValues[index];
    }

    const std::vector<T> &F() const
    {
        return m_fValues;
    }

    const std::vector<T> &X() const
    {
        return m_xValues;
    }

    void Add(T x, T f)
    {
        m_xValues.push_back(x);
        m_fValues.push_back(f);
    }

    class Iterator
    {
      public:
        Iterator(const FunctionData<T> &function, std::size_t index) : m_Function(function), m_Index(index)
        {
        }

        // Dereference operator
        std::pair<T, T> operator*() const
        {
            return {m_Function.X(m_Index), m_Function.F(m_Index)};
        }

        // Prefix increment operator
        Iterator &operator++()
        {
            ++m_Index;
            return *this;
        }

        // Postfix increment operator
        Iterator operator++(int)
        {
            Iterator temp = *this;
            ++m_Index;
            return temp;
        }

        // Equality operator
        bool operator==(const Iterator &other) const
        {
            return m_Index == other.m_Index;
        }

        // Inequality operator
        bool operator!=(const Iterator &other) const
        {
            return m_Index != other.m_Index;
        }

      private:
        const FunctionData &m_Function;
        std::size_t m_Index;
    };

    Iterator begin() const
    {
        return Iterator(*this, 0);
    }

    // End iterator
    Iterator end() const
    {
        return Iterator(*this, m_fValues.size());
    }

  private:
    std::vector<T> m_xValues;
    std::vector<T> m_fValues;
};

namespace FunctionUtils
{

template <typename T> FunctionData<T> Polynomial(std::vector<T> coefficients, const Utils::Range<T> &x)
{
    FunctionData<double> p;
    for (auto &&x : x)
    {
        double sum = 0;
        for (int j = 0; j < coefficients.size(); j++)
        {
            sum += coefficients[j] * pow(x, j);
        }
        p.Add(x, sum);
    }
    return p;
}

} // namespace FunctionUtils
