#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include "../input_videos/split/split_video/headers/image_0000.h"
#include "../input_videos/split/split_video/headers/image_0001.h"

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

            float mag = sqrtf(vx*vx + vy*vy);

            flow_mag_gray[i] =
                (uint8_t)(fminf(mag*50.0f,255.0f));

            if(mag > max_mag) max_mag = mag;
        }
    }

    //------------------------------------------------------------------
    // 2. Flow Richtung als Farbbild erzeugen
    //------------------------------------------------------------------

    for(int y=0;y<HEIGHT;y++)
    {
        for(int x=0;x<WIDTH;x++)
        {
            int i = idx(x,y);

            float vx = flow_x[i];
            float vy = flow_y[i];

            float mag = sqrtf(vx*vx + vy*vy);
            float ang = atan2f(vy,vx);

            float h = (ang + M_PI)/(2*M_PI);
            float v = mag / (max_mag + 1e-6f);

            float c = v;
            float hh = h*6.0f;

            float r,g,b;

            int sector = (int)hh;
            float f = hh - sector;

            float p = 0;
            float q = (1-f)*c;
            float t = f*c;

            switch(sector)
            {
                case 0: r=c; g=t; b=p; break;
                case 1: r=q; g=c; b=p; break;
                case 2: r=p; g=c; b=t; break;
                case 3: r=p; g=q; b=c; break;
                case 4: r=t; g=p; b=c; break;
                default:r=c; g=p; b=q; break;
            }

            flow_mag_img[i*3+0] = (uint8_t)((r)*255);
            //printf("%d",flow_mag_img[i*3+0]);
            flow_mag_img[i*3+1] = (uint8_t)((g)*255);
            //printf("%d",flow_mag_img[i*3+1]);
            flow_mag_img[i*3+2] = (uint8_t)((b)*255);
            //printf("%d",flow_mag_img[i*3+2]);
        }
    }



    //------------------------------------------------------------------
    // 3. PNG speichern (wie dein Sobel-Code)
    //------------------------------------------------------------------
    // for(int i; i <= NPIX*3; i++){
    //     printf("%d %d %d",flow_);

    // }
     if(stbi_write_png(
         "flow.png",
         WIDTH,
         HEIGHT,
         3,
         flow_mag_img,
         WIDTH*3))
     {
         printf("[INFO] flow.png saved\n");
     }
     else
     {
         printf("[WARN] Could not save flow.png\n");
     }

     printf("[INFO] Optical flow finished\n");
}
