#include "matrix.h"

int main() {
  Matrix<double> U_mat(3, 3);
  Matrix<double> b_vec(3, 1);

  U_mat(0, 1) = U_mat(0, 2) = U_mat(1, 1) = U_mat(2, 2) = 1;
  U_mat(0, 0) = 2;
  U_mat(1, 2) = -2;

  b_vec(0, 0) = 1;
  b_vec(0, 1) = -1;
  b_vec(0, 2) = 4;

  U_mat.print();
  b_vec.print();

  return 0;
}
