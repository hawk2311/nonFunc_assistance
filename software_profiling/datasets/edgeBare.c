#include image_harder_1.h
#define WIDTH 768
#define HEIGHT 512

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>


#define NPIX   (WIDTH * HEIGHT)

// Statische Arrays (liegen im DRAM bzw. im simulierten RAM)
static uint8_t gray[NPIX];
static uint8_t edges[NPIX];
static uint8_t overlay[NPIX * 3];   // RGB Overlay




int main() {

    //printf("[INFO] Starting edge detection...\n");

    //----------------------------------------------------------------------
    // 1. Bild in eindimensionales Array kopieren
    //----------------------------------------------------------------------
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            //gray[y * WIDTH + x] = image_data[y][x];
            gray[y * WIDTH + x] = image_data[y][x];
        }
    }

    //----------------------------------------------------------------------
    // 2. Sobel-Kernel definieren
    //----------------------------------------------------------------------
    const int Gx[3][3] = {
        {-1, 0, 1},
        {-2, 0, 2},
        {-1, 0, 1}
    };
    const int Gy[3][3] = {
        {-1, -2, -1},
         {0,  0,  0},
         {1,  2,  1}
    };

    //----------------------------------------------------------------------
    // 3. Sobel-Kantenberechnung
    //----------------------------------------------------------------------
    for (int y = 1; y < HEIGHT - 1; y++) {
        for (int x = 1; x < WIDTH - 1; x++) {

            int sumX = 0;
            int sumY = 0;

            for (int ky = -1; ky <= 1; ky++) {
                for (int kx = -1; kx <= 1; kx++) {
                    int pixel = gray[(y + ky) * WIDTH + (x + kx)];
                    sumX += pixel * Gx[ky + 1][kx + 1];
                    sumY += pixel * Gy[ky + 1][kx + 1];
                }
            }

            int magnitude = (int)sqrt((double)(sumX * sumX + sumY * sumY));
            if (magnitude > 255) magnitude = 255;

            edges[y * WIDTH + x] = (uint8_t)magnitude;
        }
    }

    //----------------------------------------------------------------------
    // 4. Overlay-Bild erzeugen (reines RGB-Array)
    //----------------------------------------------------------------------
    for (int i = 0; i < NPIX; i++) {
        uint8_t v = gray[i];

        overlay[3*i + 0] = v;
        overlay[3*i + 1] = v;
        overlay[3*i + 2] = v;

        if (edges[i] > 100) {
            overlay[3*i + 0] = 255; // rot
            overlay[3*i + 1] = 0;
            overlay[3*i + 2] = 0;
        }
    }

    //----------------------------------------------------------------------
    // 5. Ergebnis in export-Array kopieren (für Host)
    //----------------------------------------------------------------------
    //Printing the Output is not considered
    // for (int y = 0; y < HEIGHT; y++) {
    //     for (int x = 0; x < WIDTH; x++) {
    //         printf("%i\n", edges[y * WIDTH + x]);
    //     }
    // }

   

    //printf("[INFO] Done.\n");

    return 0;
}






