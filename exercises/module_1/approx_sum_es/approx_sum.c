#include <stdio.h>
#include <math.h>

/*
 * I Concetti da capire sono quelli della somma commutativa
 * e cosa succede se si sommano valori he non si dovrebbero sommare
 * NOTE: servono numeri belli grandi
 *NOTE: guarda a tipo 40000 che sono circa costanti
 * e' presente un bias costante rispetto al valore esatto
*/
#define PRECISION float 

#define M_PI 3.14159265358979323846

// Define in preprocessor flag -DPRECISION=float

#define TRUE 1
#define FALSE 0

#define MAX_STEPS 50000
#define SUM_RESULT M_PI* M_PI / 6

PRECISION basel(int inverse, int steps);

int main() {
    // printf("Expected value: %.16f \n", SUM_RESULT);
    // printf("Normal: %.16f, Inverted: %.16f \n", basel(FALSE, MAX_STEPS), basel(TRUE, MAX_STEPS));

    PRECISION seq_diff;
    PRECISION inv_diff;

    printf("N\tSequential\tInverted\n");
    for (int i = 1; i <= MAX_STEPS; i++) {
        seq_diff = fabs(SUM_RESULT - basel(FALSE, i));
        inv_diff = fabs(SUM_RESULT - basel(TRUE, i));

        printf("%d\t%.20f\t%.20f\n", i, seq_diff, inv_diff);
    }
}

PRECISION basel(int inverse, int steps) {
    PRECISION sum = 0;

    if (!inverse) {
        for (int i = 1; i <= steps; i++) {
            sum += (double) (1.0f / (i*i));
        }
    } else {
        for (int i = steps; i >= 1; i--) {
            sum += (double) (1.0f / (i*i));
        }
    }

    return sum;
}
