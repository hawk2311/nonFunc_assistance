
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include "image_data.h" // the image array: image_data[HEIGHT][WIDTH]

#define WIDTH 600
#define HEIGHT 600

int main() {
    int N = WIDTH * HEIGHT;

    fftw_complex *in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    // prepare input
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            int idx = y * WIDTH + x;
            in[idx][0] = image_data[y][x];
            in[idx][1] = 0.0;
        }
    }

    // perform 2D FFT
    fftw_plan plan_forward = fftw_plan_dft_2d(HEIGHT, WIDTH, in, out, FFTW_FORWARD, FFTW_ESTIMATE); //plan for complex dft
    fftw_execute(plan_forward);

    // write result in file 
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

    printf("FFT finished, output in 'fft_output.txt'.\n");
    return 0;
}

