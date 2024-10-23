#pragma once

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace tensor {
template <typename T>
class Tensor;
}

namespace tensor {

/*
 * This is a modified version of the matrix.h header provided in e-learning
 *
 * It extends matrix calculations to a more easier control over math vectors
 *
 * To create Tensors refer to the TensorBuilder class
 */
template <typename T>
class Tensor {
  /* ----------------- CREATION FUNCTIONS --------------------------------- */

 public:
  static Tensor<T> Vector(int n) { return Tensor<T>(n, 1); }

  static Tensor<T> Matrix(int n, int m) { return Tensor<T>(n, m); }

  static Tensor<T> SMatrix(int n) { return Tensor<T>(n, n); }

  static Tensor<T> Identity(int n) {
    Tensor<T> out(n, n);
    for (int i = 0; i < n; i++) out(i, i) = 1.0;
    return out;
  }

  static Tensor<T> Ones(int n) {
    Tensor<T> out(n, n);
    for (int i = 0; i < n; i++) out(i, i) = 1.0;
    return out;
  }

  static Tensor<T> One(int n) {
    Tensor<T> out(n, 1);
    for (int i = 0; i < n; i++) out(i) = 1.0;
    return out;
  }

  static Tensor<T> Random(int n) {
    Tensor<T> out(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) out(i, j) = rand() % 100;
    return out;
  }

  static Tensor<T> Random(int n, int m) {
    Tensor<T> out(n, m);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++) out(i, j) = rand() % 100;
    return out;
  }

  static Tensor<T> FromData(std::vector<T> data, int rows, int cols) {
    return Tensor<T>(std::move(data), rows, cols);
  }

  static Tensor<T> FromData(std::vector<T> data) {
    return FromData(std::move(data), data.size(), 1);
  }

 public:
  Tensor() : rows_(0), cols_(0) {}

  Tensor(Tensor<T> const &other)
      : rows_(other.rows_), cols_(other.cols_), data_(other.data_) {}

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

  /* ---------------------------------------------------------------------- */

  Tensor<T> &operator=(Tensor<T> const &in) {
    data_ = in.data_;
    rows_ = in.rows_;
    cols_ = in.cols_;
    return *this;
  }

  bool operator==(Tensor<T> const &b) {
    return rows_ == b.rows_ && cols_ == b.cols_ && data_ == b.data_;
  }

  bool &operator==(Tensor<T> const &b) const {
    return rows_ == b.rows_ && cols_ == b.cols_ && data_ == b.data_;
  }

  int Cols() { return cols_; }

  int Rows() { return rows_; }

  int Cols() const { return cols_; }

  int Rows() const { return rows_; }

  std::vector<T> &RawData() { return data_; }

 private:
  /*
   * Retrieve vector position from tensor coordinates
   * (Starting from 0,0 to nc,nr)
   */
  inline int Site(int i, int j) const { return i * cols_ + j; }

  Tensor(int rows, int cols)
      : rows_(rows), cols_(cols), data_(rows * cols, 0) {}

  Tensor(std::vector<T> data, int rows, int cols)
      : rows_(rows), cols_(cols), data_(data) {
    assert(data.size() == rows * cols);
  }

 private:
  // Values of the Tensor
  std::vector<T> data_;

  // Number of rows
  int rows_;
  // Number of columns
  int cols_;
};

/* ----------------- TENSOR IMPLEMENTATION ---------------------------------- */

template <typename T>
void Tensor<T>::SwapRows(int i, int j) {
  T *tmp = new T[cols_];

  // Move slices of the std::vector data on the correct column representation
  // in the std::vector format
  std::copy(tmp, &data_[Site(i, 0)], sizeof(T) * cols_);
  std::copy(&data_[Site(i, 0)], &cols_[Site(j, 0)], sizeof(T) * cols_);
  std::copy(&data_[Site(j, 0)], tmp, sizeof(T) * cols_);

  delete[] tmp;
}

template <typename T>
void Tensor<T>::LinearCombRows(int i, int j, T value, int final) {
  assert(i < rows_ && j < rows_);
  assert(final < rows_);

  for (int k = 0; k < cols_; k++) {
    data_[Site(final, k)] = data_[Site(i, k)] + value * data_[Site(j, k)];
  }
}

template <typename T>
T Tensor<T>::operator()(int i, int j) const {
  assert(i < rows_ && j < rows_);

  return data_[Site(i, j)];
}

template <typename T>
T &Tensor<T>::operator()(int i, int j) {
  assert(i < rows_ && j < cols_);

  return data_[Site(i, j)];
}

template <typename T>
T Tensor<T>::operator()(int i) const {
  // Tensor should be a vector
  assert(rows_ == 1 || cols_ == 1);

  // Check if vector is transposed
  if (cols_ == 1) {
    assert(i < rows_);

    return data_[Site(i, 0)];
  }

  assert(i < cols_);

  return data_[Site(0, i)];
}

template <typename T>
T &Tensor<T>::operator()(int i) {
  // Tensor should be a vector
  assert(rows_ == 1 || cols_ == 1);

  // Check if vector is transposed
  if (cols_ == 1) {
    assert(i < rows_);

    return data_[Site(i, 0)];
  }

  assert(i < cols_);

  return data_[Site(0, i)];
}

template <typename T>
Tensor<T> Tensor<T>::operator+(Tensor<T> const &b) const {
  // Tensors should have the same dimension
  assert((rows_ == b.rows_) && (cols_ == b.cols_));

  Tensor<T> out(rows_, cols_);

  for (int i = 0; i < rows_ * cols_; i++) out.data_[i] = data_[i] + b.data_[i];

  return out;
}

template <typename T>
Tensor<T> Tensor<T>::Dot(Tensor<T> const &a) const {
  // Check for compability of the tensors for a dot product
  assert((cols_ == a.rows_));

  Tensor<T> out(rows_, a.cols_);

  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < a.cols_; j++)
      for (int k = 0; k < cols_; k++) out(i, j) += data_[Site(i, k)] * a(k, j);
  return out;
}

template <typename T>
Tensor<T> Tensor<T>::Dagger() {
  Tensor<T> out(cols_, rows_);

  for (int i = 0; i < cols_; i++)
    for (int j = 0; j < rows_; j++) out(i, j) = conj(data_[Site(j, i)]);

  return out;
}

template <typename T>
T Tensor<T>::Norm() {
  T n = 0.0;

  for (int i = 0; i < rows_ * cols_; i++) {
    T h = std::abs(data_[i]);
    n += h * h;
  }

  return n;
}

template <typename T>
Tensor<T> Tensor<T>::operator-(Tensor<T> const &b) const {
  return *this + (b * (-1.0));
}

template <typename T>
Tensor<T> Tensor<T>::operator*(T const &d) const {
  Tensor<T> out(*this);

  for (int i = 0; i < rows_ * cols_; i++) out.data_[i] *= d;

  return out;
}

template <typename T>
Tensor<T> Tensor<T>::operator*(T d) {
  Tensor<T> out(*this);
  for (int i = 0; i < rows_ * cols_; i++) out.data_[i] *= d;

  return out;
}

}  // namespace tensor
