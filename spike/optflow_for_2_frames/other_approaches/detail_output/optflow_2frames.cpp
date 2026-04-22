
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include "../../../input/videos/split/split_video/headers/image_0000.h"
#include "../../../input/videos/split/split_video/headers/image_0001.h"

#define WIDTH 426
#define HEIGHT 240
#define NPIX (WIDTH*HEIGHT)

// Output arrays
static float flow_x[NPIX];
static float flow_y[NPIX];

//static uint8_t flow_mag_img[NPIX*3];
static uint8_t flow_mag_gray[NPIX];

static inline int idx(int x,int y)
{
    return y*WIDTH + x;
}

int main()
{
    //printf("[INFO] Computing optical flow...\n");

    const uint8_t* img1 = image_0000;
    const uint8_t* img2 = image_0001;

    float max_mag = 1e-6f;

    //------------------------------------------------------------------
    // 1. Optical Flow berechnen
    //------------------------------------------------------------------

    for(int y=1;y<HEIGHT-1;y++)
    {
        for(int x=1;x<WIDTH-1;x++)
        {
            int i = idx(x,y);

            float Ix = img1[idx(x+1,y)] - img1[idx(x-1,y)];

            float Iy = img1[idx(x,y+1)] - img1[idx(x,y-1)];

            float It = img2[i] - img1[i];

            float denom = Ix*Ix + Iy*Iy + 1e-4f;

            float vx = -Ix * It / denom;
            float vy = -Iy * It / denom;

            flow_x[i] = vx;
            flow_y[i] = vy;

            float mag = sqrtf(vx*vx + vy*vy);

            flow_mag_gray[i] =
                (uint8_t)(fminf(mag*50.0f,255.0f));

            if(mag > max_mag) max_mag = mag;
        }
    }

    
    for(int y = 0; y < HEIGHT; y++)
    {
        for(int x = 0; x < WIDTH; x++)
        {
            int i = idx(x,y);
            printf("%i\n", flow_mag_gray[i]); // korrekt für uint8_t
        }
        //printf("\n"); // neue Zeile = neue Bildzeile
    }
    
}
