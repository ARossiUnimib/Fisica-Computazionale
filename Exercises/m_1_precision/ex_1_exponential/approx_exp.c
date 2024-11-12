#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * Exp con semplicemente taylor
 */
double exp_approx(double x, int steps);

/*
 * Versione alternativa che utilizza taylor ma raggrupa usando il termine prece
 * dente invece di ricalcolarlo
 */

/**
 * Versione alternativa di exp_approx
 */
double exp_alt(double x, int steps);

/**
 * Calcola il fattoriale di un numero
 */
int factorial(int n);

int main(int argc, char const *argv[]) {
  double abs_error, estim_error;

  // interpreting argv[1] as steps
  if (argc < 2) {
    printf("Usage: %s <steps>\n", argv[0]);
    return 1;
  }

  int n = atoi(argv[1]);

  // printf("%f", exp_alt(1.0, MAX_STEPS));

  for (double x = 0.1; x < 1; x += 0.1) {
    abs_error = fabs(exp_approx(x, n) - exp(x));
    estim_error = (double)(powl(x, n + 1) / (double)factorial(n + 1));

    printf("%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\n", x, exp(x), exp_approx(x, n),
           abs_error, estim_error, fabs(abs_error - estim_error));
  }
}

double exp_approx(double x, int steps) {
  double sum = 0;
  for (int i = 0; i <= steps; i++) {
    sum += powl(x, i) / (double)factorial(i);
  }

  return sum;
}

// TODO: PER RELAZIONE: Si potrebbe discutere quale tra i due algoritmi ha un
// costo computazionale minore...
double exp_alt(double x, int steps) {
  double sum = 1;
  double prev = 1;
  double temp;

  for (int i = 1; i <= steps; i++) {
    temp = prev * (double)(x / i);
    sum += temp;
    prev = temp;
  }

  return sum;
}

int factorial(int n) {
  int prod = 1;

  for (int i = n; i > 0; i--) {
    prod *= i;
  }

  return prod;
}
