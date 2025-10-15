// #include <fftw3.h>
// #include <stdio.h>
// #include <stdlib.h>

// int main() {
//     int width = 600;   // Beispiel! musst du anpassen
//     int height = 600;  // Beispiel! musst du anpassen
//     int N = width * height;

//     // Komplexes FFT-Array
//     fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * height * (width/2+1));
//     double *reconstructed = (double*) fftw_malloc(sizeof(double) * N);

//     // 1. FFT-Werte aus Datei lesen
//     FILE *f = fopen("fft_output.txt", "r");
//     if (!f) { perror("fft_output.txt"); exit(1); }

//     int idx;
//     double real, imag;
//     int count = 0;
//     while (fscanf(f, "out[%d] = %lf + %lfi\n", &idx, &real, &imag) == 3) {
//         out[idx][0] = real; // Realteil
//         out[idx][1] = imag; // Imaginärteil
//         count++;
//     }
//     fclose(f);

//     printf("Eingelesen: %d FFT-Werte\n", count);

//     // 2. Inverse FFT
//     fftw_plan plan_backward = fftw_plan_dft_c2r_2d(height, width, out, reconstructed, FFTW_ESTIMATE);
//     fftw_execute(plan_backward);

//     // 3. Normierung
//     for (int i = 0; i < N; i++) {
//         reconstructed[i] /= N;
//     }

//     // 4. Ausgabe als Graustufenbild (PGM)
//     FILE *outimg = fopen("reconstructed.pgm", "wb");
//     fprintf(outimg, "P5\n%d %d\n255\n", width, height);
//     for (int i = 0; i < N; i++) {
//         double val = reconstructed[i];
//         if (val < 0) val = 0;
//         if (val > 255) val = 255;
//         unsigned char px = (unsigned char) val;
//         fwrite(&px, 1, 1, outimg);
//     }
//     fclose(outimg);

//     // Cleanup
//     fftw_destroy_plan(plan_backward);
//     fftw_free(out);
//     fftw_free(reconstructed);

//     printf("Fertig: rekonstruierte Datei 'reconstructed.pgm' erstellt.\n");
//     return 0;
// }

#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>


#define WIDTH 600
#define HEIGHT 600

int main() {
    int N = WIDTH * HEIGHT;

    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    double *reconstructed = (double*) fftw_malloc(sizeof(double) * N);

    // 1. FFT-Werte aus Datei einlesen
    FILE *fin = fopen("fft_output.txt", "r");
    if (!fin) { perror("fft_output.txt"); exit(1); }
    for (int i = 0; i < N; i++) {
        if (fscanf(fin, "%lf %lf", &out[i][0], &out[i][1]) != 2) {
            fprintf(stderr, "Fehler beim Einlesen der FFT-Daten\n");
            exit(1);
        }
    }
    fclose(fin);

    // 2. Inverse FFT (komplex → komplex)
    fftw_plan plan_backward = fftw_plan_dft_2d(HEIGHT, WIDTH, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_backward);

    // 3. Realteil extrahieren
    for (int i = 0; i < N; i++) {
        reconstructed[i] = in[i][0];  // noch NICHT durch N teilen
    }

    // 4. Skalierung auf 0–255
    double min = reconstructed[0], max = reconstructed[0];
    for (int i = 0; i < N; i++) {
        if (reconstructed[i] < min) min = reconstructed[i];
        if (reconstructed[i] > max) max = reconstructed[i];
    }

    FILE *pgm = fopen("reconstructed.pgm", "wb");
    if (!pgm) { perror("reconstructed.pgm"); exit(1); }
    fprintf(pgm, "P5\n%d %d\n255\n", WIDTH, HEIGHT);

    for (int i = 0; i < N; i++) {
        double val = (reconstructed[i] - min) / (max - min) * 255.0;
        if (val < 0) val = 0;
        if (val > 255) val = 255;
        unsigned char px = (unsigned char) val;
        fwrite(&px, 1, 1, pgm);
    }
    fclose(pgm);

    fftw_destroy_plan(plan_backward);
    fftw_free(in);
    fftw_free(out);
    fftw_free(reconstructed);

    printf("Rekonstruktion fertig: 'reconstructed.pgm' erstellt.\n");
    return 0;
}
