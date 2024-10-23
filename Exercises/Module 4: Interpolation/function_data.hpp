#pragma once

#include <cassert>
#include <cmath>
#include <cwchar>
#include <vector>

namespace func {

/**
 * @brief Class to store data of a function
 *
 * @tparam T precision
 */
template <typename T>
class FunctionData {
   public:
    FunctionData() : x_values_(), f_values_() {}

    FunctionData(std::size_t size) {
        x_values_.reserve(size);
        f_values_.reserve(size);
    }

    FunctionData(std::vector<T> x_values, std::vector<T> f_values)
        : x_values_(std::move(x_values)), f_values_(std::move(f_values)) {
        assert(x_values_.size() == f_values_.size());
    }

    std::size_t Size() const { return x_values_.size(); }

    T X(std::size_t index) const { return x_values_[index]; }

    T F(std::size_t index) const { return f_values_[index]; }

    const std::vector<T> &F() const { return f_values_; }

    const std::vector<T> &X() const { return x_values_; }

    void Add(T x, T f) {
        x_values_.push_back(x);
        f_values_.push_back(f);
    }

    class Iterator {
       public:
        Iterator(const FunctionData<T> &function, std::size_t index)
            : function_(function), index_(index) {}

        // Dereference operator
        std::pair<T, T> operator*() const {
            return {function_.X(index_), function_.F(index_)};
        }

        Iterator &operator++() {
            ++index_;
            return *this;
        }

        Iterator operator++(int) {
            Iterator temp = *this;
            ++index_;
            return temp;
        }

        bool operator==(const Iterator &other) const {
            return index_ == other.index_;
        }

        // Inequality operator
        bool operator!=(const Iterator &other) const {
            return index_ != other.index_;
        }

       private:
        const FunctionData &function_;
        std::size_t index_;
    };

    Iterator begin() const { return Iterator(*this, 0); }

    // End iterator
    Iterator end() const { return Iterator(*this, f_values_.size()); }

   private:
    std::vector<T> x_values_;
    std::vector<T> f_values_;
};

}  // namespace func
