#include <stdio.h>

int main() {

    // NOTE: la somma e' commutativa per valori in questo range di decimali
    float x = (0.7 + 0.1) + 0.3;
    float y = 0.7 + (0.1 + 0.3);

    printf("Left: %.16f, Right: %.16f \n", x, y);

    float xt = 1.e20;
    float yt = -1.e20;
    float zt = 1;

    // NOTE: la somma non e' piu' commutativa per valori cosi' alti e' necessaria
    // una maggiore precisione (per esempio una quadrupla precisione) per far 
    // ritornare il caso come alla nota precedente
    printf("Left: %.16f, Right: %.16f \n", (xt + yt) + zt, xt + (yt + zt));
}

