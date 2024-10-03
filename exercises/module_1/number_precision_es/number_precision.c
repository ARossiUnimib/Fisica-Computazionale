#include <stdio.h>

// TODO: Si discutera' il fatto che la moltiplicazione e l'addizione propagano
// diversamente: il metodo dei floating-point e' molto preciso per moltiplicazioni
// ma non per addizioni.
//
// NOTE: Inoltre si nota che valori precisi (come inf) appaiono. Essi hanno un
// valore preciso nello sstandard IEE...
int main() {
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
        printf("%.20e, %.20e \n", d, (double) (1+d));
    }

    f = 1e-13;
    printf("Single precision division: \n");

    for (int i = 0; i < 24; i++) {
        f /= 2;
        printf("%.20e, %.20e \n", f, (float) (1+f));
    }
}
