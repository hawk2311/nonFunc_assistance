// This Version calculates the average movement in x- and y- direction
//So far it is not possible to gurantee the correctness of the results


#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include "input_videos/split/split_video/headers/image_0001.h"
#include "input_videos/split/split_video/headers/image_0002.h"

#define WIDTH 1280
#define HEIGHT 720
#define NPIX (WIDTH*HEIGHT)

// Output arrays
static float flow_x[NPIX];
static float flow_y[NPIX];

static uint8_t flow_mag_img[NPIX*3];
static uint8_t flow_mag_gray[NPIX];

static inline int idx(int x,int y)
{
    return y*WIDTH + x;
}

int main()
{
    printf("[INFO] Computing optical flow...\n");

    const uint8_t* img1 = image_0001;
    const uint8_t* img2 = image_0002;

    float max_mag = 1e-6f;

    // Summen für Durchschnitt
    double sum_vx = 0.0;
    double sum_vy = 0.0;
    double sum_mag = 0.0;
    int count = 0;

    //------------------------------------------------------------------
    // 1. Optical Flow berechnen
    //------------------------------------------------------------------

    for(int y=1;y<HEIGHT-1;y++)
    {
        for(int x=1;x<WIDTH-1;x++)
        {
            int i = idx(x,y);

            float Ix =
                img1[idx(x+1,y)] -
                img1[idx(x-1,y)];

            float Iy =
                img1[idx(x,y+1)] -
                img1[idx(x,y-1)];

            float It =
                img2[i] - img1[i];

            float denom = Ix*Ix + Iy*Iy + 1e-4f;

            float vx = -Ix * It / denom;
            float vy = -Iy * It / denom;

            flow_x[i] = vx;
            flow_y[i] = vy;

            // >>> NEU: aufsummieren
            sum_vx += vx;
            sum_vy += vy;

            float mag = sqrtf(vx*vx + vy*vy);
            sum_mag += mag;

            count++;

            flow_mag_gray[i] =
                (uint8_t)(fminf(mag*50.0f,255.0f));

            if(mag > max_mag) max_mag = mag;
        }
    }

    //  Durchschnitt berechnen
    double avg_vx = sum_vx / count;
    double avg_vy = sum_vy / count;
    double avg_mag = sum_mag / count;

    printf("\n[RESULT] Average optical flow:\n");
    printf("vx = %f\n", avg_vx);
    printf("vy = %f\n", avg_vy);
    printf("magnitude = %f\n\n", avg_mag);

    
}
