// #include <stdint.h>
// #include <cmath>
// #include "images_index.h"

// #define FLOW_PAIRS (IMG_COUNT-1)

// static int16_t flow_all[FLOW_PAIRS][IMG_WIDTH*IMG_HEIGHT*2];

// static inline int idx(int x,int y)
// {
//     return y*IMG_WIDTH + x;
// }

// void compute_flow(
//     const uint8_t* img1,
//     const uint8_t* img2,
//     int16_t* flow
// )
// {
//     for(int y=1;y<IMG_HEIGHT-1;y++)
//     {
//         for(int x=1;x<IMG_WIDTH-1;x++)
//         {
//             int i = idx(x,y);

//             int Ix =
//                 img1[idx(x+1,y)] -
//                 img1[idx(x-1,y)];

//             int Iy =
//                 img1[idx(x,y+1)] -
//                 img1[idx(x,y-1)];

//             int It =
//                 img2[i] - img1[i];

//             float denom = Ix*Ix + Iy*Iy + 1e-4f;

//             float vx = -Ix * It / denom;
//             float vy = -Iy * It / denom;

//             int fi = i*2;

//             flow[fi+0] = (int16_t)(vx*32);
//             flow[fi+1] = (int16_t)(vy*32);
//         }
//     }
// }


// int main()
// {
//     for(int i=0;i<FLOW_PAIRS;i++)
//     {
//         const uint8_t* frame1 = frames[i];
//         const uint8_t* frame2 = frames[i+1];

//         compute_flow(frame1, frame2, flow_all[i]);

//         printf("Computed flow for frames %d -> %d\n", i, i+1);
//     }

//     printf("All optical flow pairs computed\n");

//     return 0;
// }

#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "input_videos/split/split_video/headers/images_index.h"

#define FLOW_PAIRS (IMG_COUNT-1)

static inline int idx(int x,int y)
{
    return y*IMG_WIDTH + x;
}

void save_ppm(const char* filename, uint8_t* img)
{
    FILE* f = fopen(filename,"wb");

    fprintf(f,"P6\n%d %d\n255\n",IMG_WIDTH,IMG_HEIGHT);
    fwrite(img,1,IMG_WIDTH*IMG_HEIGHT*3,f);

    fclose(f);
}

void compute_flow_and_visualize(
    const uint8_t* img1,
    const uint8_t* img2,
    uint8_t* vis_out
)
{
    float* vx_field = (float*)malloc(sizeof(float)*IMG_WIDTH*IMG_HEIGHT);
    float* vy_field = (float*)malloc(sizeof(float)*IMG_WIDTH*IMG_HEIGHT);

    float max_mag = 0.0f;

    // Flow berechnen
    for(int y=1;y<IMG_HEIGHT-1;y++)
    {
        for(int x=1;x<IMG_WIDTH-1;x++)
        {
            int i = idx(x,y);

            int Ix =
                img1[idx(x+1,y)] -
                img1[idx(x-1,y)];

            int Iy =
                img1[idx(x,y+1)] -
                img1[idx(x,y-1)];

            int It =
                img2[i] - img1[i];

            float denom = Ix*Ix + Iy*Iy + 1e-4f;

            float vx = -Ix * It / denom;
            float vy = -Iy * It / denom;

            vx_field[i] = vx;
            vy_field[i] = vy;

            float mag = sqrtf(vx*vx + vy*vy);
            if(mag > max_mag) max_mag = mag;
        }
    }

    // Visualisierung erzeugen
    for(int y=0;y<IMG_HEIGHT;y++)
    {
        for(int x=0;x<IMG_WIDTH;x++)
        {
            int i = idx(x,y);

            float vx = vx_field[i];
            float vy = vy_field[i];

            float mag = sqrtf(vx*vx + vy*vy);
            float ang = atan2f(vy,vx);

            // HSV → RGB (manuell)
            float h = (ang + M_PI) / (2*M_PI);
            float s = 1.0f;
            float v = mag / (max_mag + 1e-6f);

            float c = v * s;
            float xcol = c * (1 - fabs(fmod(h*6,2)-1));
            float m = v - c;

            float r,g,b;

            int sector = (int)(h*6);

            switch(sector)
            {
                case 0: r=c; g=xcol; b=0; break;
                case 1: r=xcol; g=c; b=0; break;
                case 2: r=0; g=c; b=xcol; break;
                case 3: r=0; g=xcol; b=c; break;
                case 4: r=xcol; g=0; b=c; break;
                default:r=c; g=0; b=xcol; break;
            }

            int base = i*3;

            vis_out[base+0] = (uint8_t)((r+m)*255);
            vis_out[base+1] = (uint8_t)((g+m)*255);
            vis_out[base+2] = (uint8_t)((b+m)*255);
        }
    }

    free(vx_field);
    free(vy_field);
}

int main()
{
    uint8_t* vis =
        (uint8_t*)malloc(IMG_WIDTH*IMG_HEIGHT*3);

    for(int i=0;i<FLOW_PAIRS;i++)
    {
        compute_flow_and_visualize(
            frames[i],
            frames[i+1],
            vis
        );

        char filename[128];
        sprintf(filename,
            "flow_vis_%04d_%04d.ppm",
            i,i+1);

        save_ppm(filename,vis);

        printf("Saved %s\n",filename);
    }

    free(vis);

    return 0;
}

