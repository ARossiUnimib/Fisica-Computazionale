#include <math.h>
#include <stdio.h>

/*
 * I Concetti da capire sono quelli della somma commutativa
 * e cosa succede se si sommano valori he non si dovrebbero sommare
 * NOTE: servono numeri belli grandi
 *NOTE: guarda a tipo 40000 che sono circa costanti
 * e' presente un bias costante rispetto al valore esatto
 */
#define PRECISION double

#define TRUE 1
#define FALSE 0

#define MAX_STEPS 1000000
#define SUM_RESULT \
  1.6449340668482264364724151666460251892189499012067984377355582293700074704032008738336289006197587053040043189623371906796287246870050077879351029463308662768317333093677626050952510068721400547968116

PRECISION basel(int inverse, int steps);

int main() {
  // printf("Expected value: %.16f \n", SUM_RESULT);
  // printf("Normal: %.16f, Inverted: %.16f \n", basel(FALSE, MAX_STEPS),
  // basel(TRUE, MAX_STEPS));

  PRECISION seq_diff;
  PRECISION inv_diff;

  printf("N\tSequential\tInverted\n");
  for (int i = 1; i <= MAX_STEPS; i++) {
    seq_diff = fabs(SUM_RESULT - basel(FALSE, i));
    inv_diff = fabs(SUM_RESULT - basel(TRUE, i));

    printf("%d\t%.20f\t%.20f\t%.20f\t%.20f\n", i, basel(FALSE, i),
           basel(TRUE, i), seq_diff, inv_diff);
  }
}

PRECISION basel(int inverse, int steps) {
  PRECISION sum = 0;

  if (!inverse) {
    for (int i = 1; i <= steps; i++) {
      sum += (double)(1.0f / (i * i));
    }
  } else {
    for (int i = steps; i >= 1; i--) {
      sum += (double)(1.0f / (i * i));
    }
  }

  return sum;
}
