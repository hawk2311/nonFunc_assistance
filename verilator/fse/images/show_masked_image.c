#include <stdio.h>
#include "mask_image.h"   // enthält z.B. static int mask_image[600][600]

#define SIZE 200

void write_pgm(int data[SIZE][SIZE], const char *filename) {
    FILE *f = fopen(filename, "wb");
    if (!f) {
        perror("fopen");
        return;
    }

    // PGM-Header (P5 = binäres Graustufenformat)
    fprintf(f, "P5\n%d %d\n255\n", SIZE, SIZE);

    // Pixel-Daten als Bytes schreiben
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            unsigned char pixel = (unsigned char)data[i][j];
            fwrite(&pixel, 1, 1, f);
        }
    }

    fclose(f);
}

int main(void) {
    write_pgm(mask_image, "output.pgm");
    return 0;
}
