#pragma once

#include <cmath>
#include <functional>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace Utils
{

/**
 * @brief Range struct
 *
 * @tparam T
 */
template <typename T> struct Range;

/**
 * @brief Slice a vector
 *
 * @tparam T
 * @param v
 * @param slice_range
 * @return std::vector<T>
 */
template <typename T> std::vector<T> slice(std::vector<T> const &v, Range<T> slice_range);

} // namespace Utils

namespace Utils
{

template <typename T> struct Range
{
    enum class IntervalType
    {
        Inclusive,          // Both Start and End are inclusive
        Exclusive,          // Both Start and End are exclusive
        SemiInclusiveStart, // Start inclusive, End exclusive
        SemiInclusiveEnd    // Start exclusive, End inclusive
    } IntType;

    enum class NodeType
    {
        Chebyshev,
        Fixed
    } Type;

    class Iterator
    {
        const std::vector<T> &nodes;
        size_t index;

      public:
        Iterator(const std::vector<T> &nodes, size_t index = 0) : nodes(nodes), index(index)
        {
        }

        T operator*() const
        {
            return nodes[index];
        }

        Iterator &operator++()
        {
            ++index;
            return *this;
        }

        // Postfix increment
        Iterator operator++(int)
        {
            Iterator temp = *this;
            ++(*this);
            return temp;
        }

        // Comparison operator
        bool operator!=(const Iterator &other) const
        {
            return index != other.index;
        }
    };

    T Start;
    T End;
    T Step;
    size_t NumNodes;
    std::vector<T> nodes;

    static Range Chebyshev(T start, T end, size_t numNodes, IntervalType intType = IntervalType::Inclusive)
    {
        return Range(start, end, numNodes, intType);
    }

    static Range Fixed(T start, T end, T step = 1, IntervalType intType = IntervalType::Inclusive)
    {
        return Range(start, end, step, intType);
    }

    static Range FixedNum(T start, T end, T numNodes, IntervalType intType = IntervalType::Inclusive)
    {
        return Range(start, end, ((end - start) / numNodes), intType);
    }

    Range() = default;

  private:
    // Constructor for Chebyshev nodes
    Range(T start, T end, size_t numNodes, IntervalType intType = IntervalType::Inclusive)
        : Start(start), End(end), NumNodes(numNodes), Step(0), Type(NodeType::Chebyshev), IntType(intType)
    {
        nodes.reserve(numNodes);  // Reserve space for Chebyshev nodes
        generateChebyshevNodes(); // Precompute Chebyshev nodes
    }

    // Constructor for Fixed distance nodes
    Range(T start, T end, T step, IntervalType intType = IntervalType::Inclusive)
        : Start(start), End(end), NumNodes(0), Step(step), Type(NodeType::Fixed), IntType(intType)
    {
        generateFixedDistanceNodes(); // Precompute fixed distance nodes
    }

    // Generate Chebyshev nodes
    void generateChebyshevNodes()
    {
        // Inverting the formula for Chebyshev nodes
        for (size_t k = NumNodes; k-- > 0;)
        {
            T node = (Start + End) / 2 + (End - Start) / 2 * std::cos((2 * k + 1) * M_PI / (2 * NumNodes));
            nodes.push_back(node);
        }
    }

    // Generate Fixed distance nodes
    void generateFixedDistanceNodes()
    {
        for (T current = Start; (IntType == IntervalType::Inclusive) ? (current <= End) : (current < End);
             current += Step)
        {
            nodes.push_back(current);
        }
    }

  public:
    // Begin iterator
    Iterator begin()
    {
        return Iterator(nodes);
    }

    // End iterator
    Iterator end()
    {
        return Iterator(nodes, nodes.size()); // Point to past the last element
    }

    // Const begin iterator
    Iterator begin() const
    {
        return Iterator(nodes);
    }

    // Const end iterator
    Iterator end() const
    {
        return Iterator(nodes, nodes.size()); // Point to past the last element
    }
};

template <typename T> std::vector<T> slice(std::vector<T> const &v, Range<T> slice_range)
{
    auto first = v.cbegin() + slice_range.Start;
    auto last = v.cbegin() + slice_range.End + 1;
    return std::vector<T>(first, last);
}

} // namespace Utils
