#include <stdio.h>

int main() {

    float x = (0.7 + 0.1) + 0.3;
    float y = 0.7 + (0.1 + 0.3);

    printf("Left: %.16f, Right: %.16f \n", x, y);

    float xt = 1.e20;
    float yt = -1.e20;
    float zt = 1;

    printf("Left: %.16f, Right: %.16f \n", (xt + yt) + zt, xt + (yt + zt));
}

