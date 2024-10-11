#pragma once

#include <vector>

namespace Utils
{

template <typename T> struct Range
{
    struct Iterator
    {
        T current;
        T steps;

        Iterator(T start, T step) : current(start), steps(step)
        {
        }

        // Dereference operator
        T operator*() const
        {
            return current;
        }

        // Prefix increment
        Iterator &operator++()
        {
            current += steps;
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
            return current <= other.current;
        }
    };

    T Start;
    T End;
    T Steps = 1;

    Range() = default;

    Range(T start, T end) : Start(start), End(end)
    {
    }

    Range(T start, T end, T step) : Start(start), End(end), Steps(step)
    {
    }

    // Begin iterator
    Iterator begin() const
    {
        return Iterator(Start, Steps);
    }

    // End iterator
    Iterator end() const
    {
        return Iterator(End, Steps);
    }
};

template <typename T> std::vector<T> slice(std::vector<T> const &v, Range<T> slice_range)
{
    auto first = v.cbegin() + slice_range.Start;
    auto last = v.cbegin() + slice_range.End + 1;
    return std::vector<T>(first, last);
}

} // namespace Utils
