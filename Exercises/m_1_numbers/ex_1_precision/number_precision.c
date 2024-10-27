#include <stdio.h>

int main() {
#ifdef EXERCISE_FORMAT
  // Exact exercise request

  float f = 1.2e34;
  printf("Single precision multiplication: \n");

  for (int i = 0; i < 24; i++) {
    f *= 2;
    printf("%.20e\n", f);
  }

  double d = 1.2e304;
  printf("Double precision multiplication: \n");

  for (int i = 0; i < 24; i++) {
    d *= 2;
    printf("%.20e\n", d);
  }

  d = 1e-13;
  printf("Double precision division: \n");

  for (int i = 0; i < 24; i++) {
    d /= 2;
    printf("%.20e, %.20e \n", d, (double)(1 + d));
  }

  f = 1e-13;
  printf("Single precision division: \n");

  for (int i = 0; i < 24; i++) {
    f /= 2;
    printf("%.20e, %.20e \n", f, (float)(1 + f));
  }
#else
  // Rewriting for gnuplot parsing

  float f_mult = 1.2e34;
  double d_mult = 1.2e304;
  float f_div = 1e-13;
  double d_div = 1e-13;

  printf("f_mult\td_mult\td_div\t1+d_div\tf_div\t1+f_div\n");

  for (int i = 0; i < 24; i++) {
    f_mult *= 2;
    d_mult *= 2;
    f_div /= 2;
    d_div /= 2;

    printf("%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\n", f_mult, d_mult, d_div,
           1 + d_div, f_div, 1 + f_div);
  }
#endif
}
