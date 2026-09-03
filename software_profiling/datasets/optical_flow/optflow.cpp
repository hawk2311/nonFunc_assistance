

//This code is taken from directory spike/optical_flow, there you can see the "original".

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

#include "../../input/videos/split/split_video/headers/images_index.h"

#define WIDTH IMG_WIDTH
#define HEIGHT IMG_HEIGHT
#define NPIX (WIDTH*HEIGHT)

// Flow
static float flow_x[NPIX];
static float flow_y[NPIX];



static inline int idx(int x,int y)
{
    return y*WIDTH + x;
}


int main()
{
    const uint8_t* img1;
    const uint8_t* img2;

    double sum_vx = 0.0;
    double sum_vy = 0.0;
    long count = 0;

    float threshold = 0.3f;

    //--------------------------------------------------
    // for every frame
    //--------------------------------------------------
    for(int f = 0; f < IMG_COUNT - 1; f++)
    {
        img1 = frames[f];
        img2 = frames[f+1];

        //-------------------------------
        // calculate flow
        //-------------------------------
        for(int y=1;y<HEIGHT-1;y++)
        {
            for(int x=1;x<WIDTH-1;x++)
            {
                int i = idx(x,y);

                float Ix = img1[idx(x+1,y)] - img1[idx(x-1,y)];
                float Iy = img1[idx(x,y+1)] - img1[idx(x,y-1)];
                float It = img2[i] - img1[i];

                float denom = Ix*Ix + Iy*Iy + 1e-4f;

                flow_x[i] = -Ix * It / denom;
                flow_y[i] = -Iy * It / denom;
            }
        }

        //-------------------------------
        // gather values for average
        //-------------------------------
        for(int i = 0; i < NPIX; i++)
        {
            float vx = flow_x[i];
            float vy = flow_y[i];

            float mag = sqrtf(vx*vx + vy*vy);

            if(mag > threshold)
            {
                sum_vx += vx;
                sum_vy += vy;
                count++;
            }
        }
    }

    //--------------------------------------------------
    // calculate average 
    //--------------------------------------------------
    float avg_vx = 0.0f;
    float avg_vy = 0.0f;

    if(count > 0)
    {
        avg_vx = sum_vx / count;
        avg_vy = sum_vy / count;
    }

    //printf("Average Flow: vx=%f vy=%f\n", avg_vx, avg_vy);

   

    

    

    return 0;
}


