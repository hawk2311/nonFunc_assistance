//Use of kissFFT --> working, executable with fft_test

// #include <stdio.h>
// #include <stdint.h>
// #include "kissfft/kiss_fft.h"
// #include "image_data.h"


// #define N 64 //FFT length with 64 values

// int main() {
//     kiss_fft_cfg cfg = kiss_fft_alloc(N, 0, NULL, NULL); // 0 = Forward FFT
//     kiss_fft_cpx in[N], out[N];

//     // Erste Zeile des Bildes in den FFT-Input kopieren
//     for (int i = 0; i < N; i++) {
//         in[i].r = image_data[0][i]; // nur reale Anteile
//         in[i].i = 0;
//     }

//     // FFT ausführen
//     kiss_fft(cfg, in, out);

//     // Ergebnis ausgeben
//     for (int i = 0; i < N; i++) {
//         printf("out[%2d] = %.2f + %.2fi\n", i, out[i].r, out[i].i);
//     }

//     free(cfg);
//     return 0;
// }

//------------------------------------------------------------------------------------------------------

//Use of FFTW --> executeable with fftw

// #include <stdio.h>
// #include <stdint.h>
// #include <fftw3.h> 
// #include "image_data.h"

// #define N 600 //it was 64 before

// int main() {
//     // Eingabe- und Ausgabe-Arrays
//     fftw_complex in[N], out[N];

//     // Plan für Vorwärts-FFT erstellen
//     fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

//     // Erste Zeile des Bildes in den FFT-Input kopieren
//     for (int i = 0; i < N; i++) {
//         in[i][0] = image_data[0][i]; // Reale Anteile
//         in[i][1] = 0.0;              // Imaginäre Anteile auf 0
//     }

//     // FFT ausführen
//     fftw_execute(plan);

//     // Ergebnis ausgeben
//     for (int i = 0; i < N; i++) {
//         printf("out[%2d] = %.2f + %.2fi\n", i, out[i][0], out[i][1]);
//     }

//     // Speicher freigeben
//     fftw_destroy_plan(plan);
//     fftw_cleanup();

//     return 0;
// }

//--------------------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include "image_data.h" // Dein Bildarray: image_data[HEIGHT][WIDTH]

#define WIDTH 600
#define HEIGHT 600

int main() {
    int N = WIDTH * HEIGHT;

    fftw_complex *in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    // Input vorbereiten
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            int idx = y * WIDTH + x;
            in[idx][0] = image_data[y][x];
            in[idx][1] = 0.0;
        }
    }

    // 2D FFT
    fftw_plan plan_forward = fftw_plan_dft_2d(HEIGHT, WIDTH, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    // FFT-Ergebnis in Datei schreiben
    FILE *fout = fopen("fft_output.txt", "w");
    if (!fout) { perror("fft_output.txt"); exit(1); }
    for (int i = 0; i < N; i++) {
        fprintf(fout, "%.10f %.10f\n", out[i][0], out[i][1]);
    }
    fclose(fout);

    //     for (int i = 0; i < N; i++) {
    //     printf("out[%2d] = %.2f + %.2fi\n", i, out[i][0], out[i][1]);
    // }

    fftw_destroy_plan(plan_forward);
    fftw_free(in);
    fftw_free(out);

    printf("FFT fertig, Ausgabe in 'fft_output.txt'.\n");
    return 0;
}

