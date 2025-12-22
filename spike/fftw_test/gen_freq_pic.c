#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#define WIDTH 600
#define HEIGHT 600

int main() {
    int N = WIDTH * HEIGHT;
    double *magnitude = malloc(sizeof(double) * N);
    double re, im, max_val = 0.0;

    FILE *fin = fopen("fft_output.txt", "r");
    if (!fin) { perror("fft_output.txt"); exit(1); }

    // FFT-Ergebnis aus Datei einlesen und Magnitude berechnen
    for (int i = 0; i < N; i++) {
        if (fscanf(fin, "%lf %lf", &re, &im) != 2) {
            fprintf(stderr, "Error while reading the line %d\n", i);
            exit(1);
        }
        magnitude[i] = log(1.0 + sqrt(re * re + im * im));
        if (magnitude[i] > max_val)
            max_val = magnitude[i];
    }
    fclose(fin);

    // Nullfrequenz ins Zentrum verschieben (fftshift)
    double *shifted = malloc(sizeof(double) * N);
    int halfW = WIDTH / 2, halfH = HEIGHT / 2;
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            int srcY = (y + halfH) % HEIGHT;
            int srcX = (x + halfW) % WIDTH;
            shifted[y * WIDTH + x] = magnitude[srcY * WIDTH + srcX];
        }
    }

    // In 0–255 skalieren und als PGM speichern
    FILE *fout = fopen("fft_visualized.pgm", "w");
    if (!fout) { perror("fft_visualized.pgm"); exit(1); }

    fprintf(fout, "P2\n%d %d\n255\n", WIDTH, HEIGHT);
    for (int i = 0; i < N; i++) {
        int val = (int)(255.0 * shifted[i] / max_val);
        fprintf(fout, "%d\n", val);
    }
    fclose(fout);

    free(magnitude);
    free(shifted);

    printf("Frequenzbild saved: fft_visualized.pgm\n");
    return 0;
}
