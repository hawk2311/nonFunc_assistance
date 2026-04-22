

// // #include <stdio.h>
// // #include <stdint.h>
// // #include <math.h>

// // #include "input_videos/split/split_video/headers/images_index.h"

// // #define WIDTH 1280
// // #define HEIGHT 720
// // #define NPIX (WIDTH*HEIGHT)
// // #define FLOW_PAIRS (IMG_COUNT-1)

// // // Output arrays
// // static float flow_x[NPIX];
// // static float flow_y[NPIX];

// // static uint8_t flow_mag_img[NPIX*3];
// // static uint8_t flow_mag_gray[NPIX];

// // static inline int idx(int x,int y)
// // {
// //     return y*WIDTH + x;
// // }

// // void save_ppm(const char* filename, uint8_t* img)
// // {
// //     FILE* f = fopen(filename,"wb");

// //     if(!f)
// //     {
// //         printf("[WARN] Could not open %s\n",filename);
// //         return;
// //     }

// //     fprintf(f,"P6\n%d %d\n255\n",WIDTH,HEIGHT);

// //     fwrite(img,1,NPIX*3,f);

// //     fclose(f);

// //     printf("[INFO] %s saved\n",filename);
// // }

// // void calc_optflow(const uint8_t* img1, const uint8_t* img2, int id)
// // {
// //     float max_mag = 1e-6f;

// //     for(int y=1;y<HEIGHT-1;y++)
// //     {
// //         for(int x=1;x<WIDTH-1;x++)
// //         {
// //             int i = idx(x,y);

// //             float Ix = img1[idx(x+1,y)] - img1[idx(x-1,y)];
// //             float Iy = img1[idx(x,y+1)] - img1[idx(x,y-1)];
// //             float It = img2[i] - img1[i];

// //             float denom = Ix*Ix + Iy*Iy + 1e-4f;

// //             float vx = -Ix * It / denom;
// //             float vy = -Iy * It / denom;

// //             flow_x[i] = vx;
// //             flow_y[i] = vy;

// //             float mag = sqrtf(vx*vx + vy*vy);

// //             flow_mag_gray[i] =
// //                 (uint8_t)(fminf(mag*50.0f,255.0f));

// //             if(mag > max_mag) max_mag = mag;
// //         }
// //     }

// //     for(int y=0;y<HEIGHT;y++)
// //     {
// //         for(int x=0;x<WIDTH;x++)
// //         {
// //             int i = idx(x,y);

// //             float vx = flow_x[i];
// //             float vy = flow_y[i];

// //             float mag = sqrtf(vx*vx + vy*vy);
// //             float ang = atan2f(vy,vx);

// //             float h = (ang + M_PI)/(2*M_PI);
// //             float v = mag/(max_mag + 1e-6f);

// //             float c = v;
// //             float hh = h*6.0f;

// //             float r,g,b;

// //             int sector = (int)hh;
// //             float f = hh - sector;

// //             float p = 0;
// //             float q = (1-f)*c;
// //             float t = f*c;

// //             switch(sector)
// //             {
// //                 case 0: r=c; g=t; b=p; break;
// //                 case 1: r=q; g=c; b=p; break;
// //                 case 2: r=p; g=c; b=t; break;
// //                 case 3: r=p; g=q; b=c; break;
// //                 case 4: r=t; g=p; b=c; break;
// //                 default:r=c; g=p; b=q; break;
// //             }

// //             flow_mag_img[i*3+0] = (uint8_t)(r*255);
// //             flow_mag_img[i*3+1] = (uint8_t)(g*255);
// //             flow_mag_img[i*3+2] = (uint8_t)(b*255);
// //         }
// //     }

// //     char filename[64];
// //     sprintf(filename,"flow_%04d.ppm",id);

// //     save_ppm(filename,flow_mag_img);

// //     printf("[INFO] Optical flow finished for frame %d\n",id);
// // }

// // int main()
// // {
// //     printf("[INFO] Computing optical flow...\n");

// //     for(int i=0;i<FLOW_PAIRS;i++)
// //     {
// //         const uint8_t* frame1 = frames[i];
// //         const uint8_t* frame2 = frames[i+1];

// //         calc_optflow(frame1,frame2,i);

// //         printf("Computed flow for frames %d -> %d\n",i,i+1);
// //     }

// //     printf("All optical flow pairs computed\n");

// //     return 0;
// // }


// #include <stdio.h>
// #include <stdint.h>
// #include <math.h>

// #include "input_videos/split/split_video/headers/images_index.h"

// #define WIDTH IMG_WIDTH
// #define HEIGHT IMG_HEIGHT
// #define NPIX (WIDTH*HEIGHT)

// // Output arrays
// static float flow_x[NPIX];
// static float flow_y[NPIX];

// static inline int idx(int x,int y)
// {
//     return y*WIDTH + x;
// }

// int main()
// {
//     printf("[INFO] Computing optical flow for %d frames...\n\n", IMG_COUNT);

//     // >>> Global über gesamtes Video
//     double global_sum_vx = 0.0;
//     double global_sum_vy = 0.0;
//     double global_sum_mag = 0.0;
//     int global_count = 0;

//     //------------------------------------------------------------------
//     // Loop über alle Frame-Paare
//     //------------------------------------------------------------------
//     for(int f = 0; f < IMG_COUNT - 1; f++)
//     {
//         const uint8_t* img1 = frames[f];
//         const uint8_t* img2 = frames[f+1];

//         double sum_vx = 0.0;
//         double sum_vy = 0.0;
//         double sum_mag = 0.0;
//         int count = 0;

//         //------------------------------------------------------------------
//         // Optical Flow für dieses Frame-Paar
//         //------------------------------------------------------------------
//         for(int y=1;y<HEIGHT-1;y++)
//         {
//             for(int x=1;x<WIDTH-1;x++)
//             {
//                 int i = idx(x,y);

//                 float Ix =
//                     img1[idx(x+1,y)] -
//                     img1[idx(x-1,y)];

//                 float Iy =
//                     img1[idx(x,y+1)] -
//                     img1[idx(x,y-1)];

//                 float It =
//                     img2[i] - img1[i];

//                 float denom = Ix*Ix + Iy*Iy + 1e-4f;

//                 float vx = -Ix * It / denom;
//                 float vy = -Iy * It / denom;

//                 flow_x[i] = vx;
//                 flow_y[i] = vy;

//                 float mag = sqrtf(vx*vx + vy*vy);

//                 // Optional: nur sinnvolle Pixel verwenden
//                 // if(Ix*Ix + Iy*Iy > 50.0f)
//                 // {
//                     sum_vx += vx;
//                     sum_vy += vy;
//                     sum_mag += mag;
//                     count++;

//                     // >>> auch global aufsummieren
//                     global_sum_vx += vx;
//                     global_sum_vy += vy;
//                     global_sum_mag += mag;
//                     global_count++;
//                 //}
//             }
//         }

//         //------------------------------------------------------------------
//         // Ergebnis für dieses Frame-Paar
//         //------------------------------------------------------------------
//         if(count > 0)
//         {
//             double avg_vx = sum_vx / count;
//             double avg_vy = sum_vy / count;
//             double avg_mag = sum_mag / count;

//             printf("[PAIR %02d -> %02d] vx = %8.4f | vy = %8.4f | mag = %8.4f\n",
//                    f, f+1, avg_vx, avg_vy, avg_mag);
//         }
//         else
//         {
//             printf("[PAIR %02d -> %02d] No valid pixels\n", f, f+1);
//         }
//     }

//     //------------------------------------------------------------------
//     // Globales Ergebnis über gesamtes Video
//     //------------------------------------------------------------------
//     printf("\n----------------------------------------\n");

//     if(global_count > 0)
//     {
//         double global_avg_vx = global_sum_vx / global_count;
//         double global_avg_vy = global_sum_vy / global_count;
//         double global_avg_mag = global_sum_mag / global_count;

//         printf("[GLOBAL RESULT]\n");
//         printf("vx = %f\n", global_avg_vx);
//         printf("vy = %f\n", global_avg_vy);
//         printf("magnitude = %f\n", global_avg_mag);
//     }
//     else
//     {
//         printf("[GLOBAL RESULT] No valid data\n");
//     }

//     printf("\n[INFO] Finished\n");

//     return 0;
// }

#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include "input_videos/split/split_video/headers/images_index.h"

#define WIDTH IMG_WIDTH
#define HEIGHT IMG_HEIGHT
#define NPIX (WIDTH*HEIGHT)
#define NPAIRS (IMG_COUNT - 1)

// >>> NEU: komplettes Flow-Array für alle Frame-Paare
static float flow_x[NPAIRS][NPIX];
static float flow_y[NPAIRS][NPIX];

static inline int idx(int x,int y)
{
    return y*WIDTH + x;
}

int main()
{
    //printf("[INFO] Computing optical flow for %d frame pairs...\n\n", NPAIRS);

    //------------------------------------------------------------------
    // Loop über alle Frame-Paare
    //------------------------------------------------------------------
    for(int f = 0; f < NPAIRS; f++)
    {
        const uint8_t* img1 = frames[f];
        const uint8_t* img2 = frames[f+1];

        //------------------------------------------------------------------
        // Optical Flow für dieses Frame-Paar
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

                // >>> NEU: Speicherung pro Frame-Paar
                flow_x[f][i] = vx;
                flow_y[f][i] = vy;
            }
        }

        //printf("[INFO] Frame pair %d -> %d processed\n", f, f+1);
    }

    //------------------------------------------------------------------
    // Gesamtes Array ausgeben
    //------------------------------------------------------------------
    //printf("\n[OUTPUT] Full optical flow data:\n");

    for(int f = 0; f < NPAIRS; f++)
    {
        //printf("\n=== Frame Pair %d -> %d ===\n", f, f+1);

        for(int i = 0; i < NPIX; i++)
        {
            printf("(%f, %f)\n", flow_x[f][i], flow_y[f][i]);
        }
    }

    //printf("\n[INFO] Finished\n");

    return 0;
}