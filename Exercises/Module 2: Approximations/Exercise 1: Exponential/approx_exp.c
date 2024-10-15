/*
 * TODO:: plottare i diversi dati in funzione di N (diverso colore) e X
 * Spiegazione del risultato dell errore che aumenta: il "leading scaling" e'
 * quello dell esercizio ma mano a mano che x non e' piu' x << 1 non e' piu'
 * considerabile come x^(n+1) ...
 */

#include <math.h>
#include <stdio.h>

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

#define MAX_STEPS 4

int main()
{
    double abs_error, estim_error;

    // printf("%f", exp_alt(1.0, MAX_STEPS));

    for (double x = 0.1; x < 1; x += 0.1)
    {
        for (int n = 0; n < MAX_STEPS; n++)
        {
            abs_error = fabs(exp_approx(x, n) - exp(x));
            estim_error = (double)(powl(x, n + 1) / (double)factorial(n + 1));

            printf("%.7f\t%.7f\t%.7f\t%.7f\n", x, abs_error, estim_error, fabs(abs_error - estim_error));
        }
    }
}

double exp_approx(double x, int steps)
{
    double sum = 0;
    for (int i = 0; i <= steps; i++)
    {
        sum += powl(x, i) / (double)factorial(i);
    }

    return sum;
}

// TODO: PER RELAZIONE: Si potrebbe discutere quale tra i due algoritmi ha un costo computazionale
// minore...
double exp_alt(double x, int steps)
{
    double sum = 1;
    double prev = 1;
    double temp;

    for (int i = 1; i <= steps; i++)
    {
        temp = prev * (double)(x / i);
        sum += temp;
        prev = temp;
    }

    return sum;
}

int factorial(int n)
{
    int prod = 1;

    for (int i = n; i > 0; i--)
    {
        prod *= i;
    }

    return prod;
}
