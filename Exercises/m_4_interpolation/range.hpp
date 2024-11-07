#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <vector>

namespace func {

template <typename T>
struct Range {
  enum class IntervalType {
    kInclusive,           // [a, b]
    kExclusive,           // (a, b)
    kSemiInclusiveStart,  // [a, b]
    kSemiInclusiveEnd     // (a, b]
  };

  enum class NodeType { kChebyshev, kFixed };

  static Range Chebyshev(T start, T end, size_t num_nodes,
                         IntervalType intrv_type = IntervalType::kInclusive) {
    return Range(start, end, num_nodes, intrv_type);
  }

  static Range Fixed(T start, T end, T step = 1,
                     IntervalType intrv_type = IntervalType::kInclusive) {
    return Range(start, end, step, intrv_type);
  }

  static Range FixedNum(T start, T end, T num_nodes,
                        IntervalType intrv_type = IntervalType::kInclusive) {
    return Range(start, end, ((end - start) / num_nodes), intrv_type);
  }

 public:
  Range() = default;

 private:
  // Constructor for Chebyshev nodes
  Range(T start, T end, size_t num_nodes,
        IntervalType intrv_type = IntervalType::kInclusive)
      : start_(start),
        end_(end),
        nodes_size_(num_nodes),
        step_(0),
        node_type_(NodeType::kChebyshev),
        intrv_type_(intrv_type) {
    nodes_.reserve(num_nodes);  // Reserve space for Chebyshev nodes
    assert(num_nodes > 0);

    GenerateChebyshevNodes();  // Precompute Chebyshev nodes
  }

  // Constructor for Fixed distance nodes
  Range(T start, T end, T step,
        IntervalType intrv_type = IntervalType::kInclusive)
      : start_(start),
        end_(end),
        nodes_size_(0),
        step_(step),
        node_type_(NodeType::kFixed),
        intrv_type_(intrv_type) {
    assert(step_ > 0);

    GenerateFixedNodes();  // Precompute fixed distance nodes
  }

  // Generate Chebyshev nodes
  void GenerateChebyshevNodes();

  // Generate Fixed distance nodes
  void GenerateFixedNodes();

  /* --------------------- MEMBERS ---------------------------------------- */

 private:
  T start_;
  T end_;

  size_t nodes_size_;
  std::vector<T> nodes_;

  IntervalType intrv_type_;
  NodeType node_type_;

  T step_;

  /* --------------------- GETTERS ---------------------------------------- */

 public:
  T Start() const { return start_; }

  T End() const { return end_; }

  size_t NodesSize() { return nodes_size_; }

  std::vector<T> &Nodes() { return nodes_; }

  NodeType NodesType() const { return node_type_; }

  T Step() const { return step_; }

  /* ------------------ ITERATOR AND EOF ---------------------------------- */

 public:
  class Iterator {
    const std::vector<T> &nodes;
    size_t index;

   public:
    Iterator(const std::vector<T> &nodes, size_t index = 0)
        : nodes(nodes), index(index) {}

    T operator*() const { return nodes[index]; }

    Iterator &operator++() {
      ++index;
      return *this;
    }

    // Postfix increment
    Iterator operator++(int) {
      Iterator temp = *this;
      ++(*this);
      return temp;
    }

    // Comparison operator
    bool operator!=(const Iterator &other) const {
      return index != other.index;
    }
  };

  // Begin iterator
  Iterator begin() { return Iterator(nodes_); }

  // End iterator
  Iterator end() {
    return Iterator(nodes_,
                    nodes_.size());  // Point to past the last element
  }

  // Const begin iterator
  Iterator begin() const { return Iterator(nodes_); }

  // Const end iterator
  Iterator end() const {
    return Iterator(nodes_,
                    nodes_.size());  // Point to past the last element
  }

  /* ---------------------------------------------------------------------- */
};

/* --------------- RANGE IMPLEMENTATION ------------------------------------- */

template <typename T>
void Range<T>::GenerateChebyshevNodes() {
  // Inverting the formula for Chebyshev nodes
  for (size_t k = nodes_size_; k-- > 0;) {
    T node =
        (start_ + end_) / 2 +
        (end_ - start_) / 2 * std::cos((2 * k + 1) * M_PI / (2 * nodes_size_));
    nodes_.push_back(node);
  }
}

template <typename T>
void Range<T>::GenerateFixedNodes() {
// Used to check which exercises are using this behaviour to avoid possible
// bugs
#ifndef NDEBUG
  if (sizeof(T) == sizeof(float) || sizeof(T) == sizeof(double)) {
    std::cout << "WARNING: using floating point can cause the last digit of "
                 "Range to be dropped"
              << std::endl;
  }
#endif

  // FIX: cause of the error propagation of floating points
  // it is not guaranteed that "current" at final step is equal to
  // "end_"!
  double _precision_fix;

  switch (sizeof(T)) {
    case sizeof(float):
      _precision_fix = 1e-5;
      break;
    case sizeof(double):
      _precision_fix = 1e-15;
      break;
    default:
      _precision_fix = 0;
      break;
  }

  for (T current = start_; (intrv_type_ == IntervalType::kInclusive)
                               ? (current <= end_ + _precision_fix)
                               : (current < end_);
       current += step_) {
    nodes_.push_back(current);
  }
}

}  // namespace func
