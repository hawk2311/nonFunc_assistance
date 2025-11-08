// #include <stdio.h>
// #include <stdint.h>
// #include <math.h>

// #define WIDTH  600
// #define HEIGHT 600
// #define FRAME_COUNT 2   // Anzahl der Frames (später anpassen)

// #include "outdir/headers/image_0000.h"
// #include "outdir/headers/image_0001.h"
// // Wenn du mehr Frames hast, einfach #include "frame2.h" usw.

// static void compute_optical_flow(const uint8_t prev[HEIGHT][WIDTH],
//                                  const uint8_t next[HEIGHT][WIDTH]) {
//     // Fenstergröße für Lucas-Kanade
//     const int W = 2;

//     printf("START_FLOW\n");

//     for (int y = W; y < HEIGHT - W; y++) {
//         for (int x = W; x < WIDTH - W; x++) {
//             float Ix = 0.0f, Iy = 0.0f, It = 0.0f;

//             // lokale Gradienten berechnen
//             for (int dy = -W; dy <= W; dy++) {
//                 for (int dx = -W; dx <= W; dx++) {
//                     float p1 = (float)prev[y + dy][x + dx];
//                     float p2 = (float)next[y + dy][x + dx];
//                     float gx = (float)prev[y + dy][x + dx + 1] - (float)prev[y + dy][x + dx - 1];
//                     float gy = (float)prev[y + dy + 1][x + dx] - (float)prev[y + dy - 1][x + dx];
//                     Ix += gx * gx;
//                     Iy += gy * gy;
//                     It += gx * (p2 - p1) + gy * (p2 - p1);
//                 }
//             }

//             // Vermeide Division durch Null
//             float denom = (Ix + Iy);
//             float u = 0.0f;
//             float v = 0.0f;

//             if (denom != 0.0f) {
//                 u = -It * Ix / denom;
//                 v = -It * Iy / denom;
//             }

//             // Ausgabe als einfache Zahlen
//             printf("%f %f\n", u, v);
//         }
//     }

//     printf("END_FLOW\n");
// }

// int main(void) {
//     printf("Lucas–Kanade Optical Flow Start\n");

//     // Beispiel: Nur Frame 0 → Frame 1
//     compute_optical_flow(image_0000, image_0001);

//     printf("Optical Flow Done\n");
//     return 0;
// }


#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "outdir/headers/images_index.h"   // Enthält alle Frames & Metadaten

#define WIN_RADIUS 2        // Lucas–Kanade Fensterhalbweite
#define EPSILON    1e-5f    // Numerische Stabilität

//-----------------------------------------------------------
// Hilfsfunktion: Optical Flow zwischen zwei Frames berechnen
//-----------------------------------------------------------
static void compute_optical_flow(const uint8_t *prev, const uint8_t *next,
                                 int width, int height, int frame_index)
{
    printf("START_FLOW %d -> %d\n", frame_index, frame_index + 1);

    for (int y = WIN_RADIUS; y < height - WIN_RADIUS; y++) {
        for (int x = WIN_RADIUS; x < width - WIN_RADIUS; x++) {
            float sum_ix2 = 0.0f, sum_iy2 = 0.0f;
            float sum_ixiy = 0.0f, sum_ixt = 0.0f, sum_iyt = 0.0f;

            for (int dy = -WIN_RADIUS; dy <= WIN_RADIUS; dy++) {
                for (int dx = -WIN_RADIUS; dx <= WIN_RADIUS; dx++) {
                    int xx = x + dx;
                    int yy = y + dy;

                    float ix = (float)prev[yy * width + (xx + 1)] - (float)prev[yy * width + (xx - 1)];
                    float iy = (float)prev[(yy + 1) * width + xx] - (float)prev[(yy - 1) * width + xx];
                    float it = (float)next[yy * width + xx] - (float)prev[yy * width + xx];

                    sum_ix2  += ix * ix;
                    sum_iy2  += iy * iy;
                    sum_ixiy += ix * iy;
                    sum_ixt  += ix * it;
                    sum_iyt  += iy * it;
                }
            }

            float det = sum_ix2 * sum_iy2 - sum_ixiy * sum_ixiy;
            float u = 0.0f, v = 0.0f;
            if (fabsf(det) > EPSILON) {
                u = (sum_iy2 * (-sum_ixt) - sum_ixiy * (-sum_iyt)) / det;
                v = (sum_ix2 * (-sum_iyt) - sum_ixiy * (-sum_ixt)) / det;
            }

            // Ausgabe pro Pixel (kann später per Python weiterverarbeitet werden)
            printf("%f %f\n", u, v);
        }
    }

    printf("END_FLOW %d -> %d\n", frame_index, frame_index + 1);
}

//-----------------------------------------------------------
// Hauptprogramm
//-----------------------------------------------------------
int main(void)
{
    printf("Lucas–Kanade Optical Flow Start (RISC-V Bare Metal)\n");
    printf("Frames: %d, Size: %dx%d\n", IMG_COUNT, IMG_WIDTH, IMG_HEIGHT);

    for (int i = 0; i < IMG_COUNT - 1; i++) {
        compute_optical_flow(frames[i], frames[i + 1], IMG_WIDTH, IMG_HEIGHT, i);
    }

    printf("Optical Flow Done\n");
    return 0;
}

