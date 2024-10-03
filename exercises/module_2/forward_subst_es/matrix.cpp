#include "matrix.h"

template <typename T> void Matrix<T>::swap_rows(int i, int j) {
  T *tmp = new T[nc];
  memcpy(tmp, &data[this->site(i, 0)], sizeof(T) * nc);
  memcpy(&data[this->site(i, 0)], &data[this->site(j, 0)], sizeof(T) * nc);
  memcpy(&data[this->site(j, 0)], tmp, sizeof(T) * nc);
}

template <typename T> void Matrix<T>::ones() {
  for (int i = 0; i < nr * nc; i++)
    data[i] = 1.0;
}

template <typename T> void Matrix<T>::random() {
  for (int i = 0; i < nr * nc; i++)
    data[i] = (T)(std::rand() / (double)RAND_MAX);
}

template <typename T> Matrix<T> &Matrix<T>::operator=(const Matrix<T> &in) {
  nr = in.nr;
  nc = in.nc;
  data = in.data;
  return *this;
}

template <typename T> Matrix<T> Matrix<T>::operator+(const Matrix<T> &b) const {
  Matrix<T> out(nr, nc);
  assert((nr == b.nr) && (nc == b.nc));
  for (int i = 0; i < nr * nc; i++)
    out.data[i] = data[i] + b.data[i];
  return out;
}

template <typename T> Matrix<T> Matrix<T>::dot(const Matrix<T> &a) const {
  Matrix<T> out(nr, a.nc);
  assert((nc == a.nr));
  for (int i = 0; i < nr; i++)
    for (int j = 0; j < a.nc; j++)
      for (int k = 0; k < nc; k++)
        out(i, j) += data[this->site(i, k)] * a(k, j);
  return out;
}

template <typename T>
template <typename T2>
Matrix<T> Matrix<T>::operator*(T2 d) const {
  Matrix<T> out(*this);
  for (int i = 0; i < nr * nc; i++)
    out.data[i] *= d;
  return out;
}

template <typename T> Matrix<T> Matrix<T>::dagger() {
  Matrix<T> out(nc, nr);
  for (int i = 0; i < nc; i++)
    for (int j = 0; j < nr; j++)
      out(i, j) = conj(data[this->site(j, i)]);
  return out;
}

template <typename T> double Matrix<T>::norm() {
  double n = 0.0;
  for (int i = 0; i < nr * nc; i++) {
    double h = std::abs(this->data[i]);
    n += h * h;
  }
  return n;
}

template <typename T> void Matrix<T>::print() {
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++)
      std::cout << data[this->site(i, j)];
    std::cout << "\n";
  }
}
