//Calculation for the whole video and drawing lines as output
//Result can be seen in flow_lines_6.png


// #include <stdio.h>
// #include <stdint.h>
// #include <math.h>
// #include <stdlib.h>

// #include "../../../input/videos/split/split_video/headers/images_index.h"

// #define WIDTH IMG_WIDTH
// #define HEIGHT IMG_HEIGHT
// #define NPIX (WIDTH*HEIGHT)

// // Flow
// static float flow_x[NPIX];
// static float flow_y[NPIX];

// // Trajektorien (merken Positionen!)
// #define MAX_POINTS ((WIDTH/15)*(HEIGHT/15))
// static float px[MAX_POINTS];
// static float py[MAX_POINTS];

// // Visualisierung
// static uint8_t flow_vis[NPIX*3];

// static inline int idx(int x,int y)
// {
//     return y*WIDTH + x;
// }

// //--------------------------------------------------
// void draw_line(uint8_t* img, int x0, int y0, int x1, int y1)
// {
//     int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
//     int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
//     int err = dx + dy;

//     while(1)
//     {
//         if(x0 >= 0 && x0 < WIDTH && y0 >= 0 && y0 < HEIGHT)
//         {
//             int i = (y0 * WIDTH + x0) * 3;
//             img[i+0] = 0;
//             img[i+1] = 255;
//             img[i+2] = 255;
//         }

//         if(x0 == x1 && y0 == y1) break;

//         int e2 = 2 * err;
//         if(e2 >= dy) { err += dy; x0 += sx; }
//         if(e2 <= dx) { err += dx; y0 += sy; }
//     }
// }

// int main()
// {
//     const uint8_t* img1;
//     const uint8_t* img2;

//     int step = 15;
//     float scale = 3.0f;

//     //------------------------------------------------------------------
//     // Hintergrund = letztes Frame
//     //------------------------------------------------------------------
//     const uint8_t* last = frames[IMG_COUNT-1];

//     for(int i = 0; i < NPIX; i++)
//     {
//         flow_vis[i*3+0] = last[i];
//         flow_vis[i*3+1] = last[i];
//         flow_vis[i*3+2] = last[i];
//     }

//     //------------------------------------------------------------------
//     // 1. Sample-Punkte initialisieren
//     //------------------------------------------------------------------
//     int pcount = 0;

//     for(int y = 0; y < HEIGHT; y += step)
//     {
//         for(int x = 0; x < WIDTH; x += step)
//         {
//             px[pcount] = x;
//             py[pcount] = y;
//             pcount++;
//         }
//     }

//     //------------------------------------------------------------------
//     // 2. Über alle Frames laufen
//     //------------------------------------------------------------------
//     for(int f = 0; f < IMG_COUNT - 1; f++)
//     {
//         img1 = frames[f];
//         img2 = frames[f+1];

//         //-------------------------------
//         // Flow EINMAL berechnen
//         //-------------------------------
//         for(int y=1;y<HEIGHT-1;y++)
//         {
//             for(int x=1;x<WIDTH-1;x++)
//             {
//                 int i = idx(x,y);

//                 float Ix = img1[idx(x+1,y)] - img1[idx(x-1,y)];
//                 float Iy = img1[idx(x,y+1)] - img1[idx(x,y-1)];
//                 float It = img2[i] - img1[i];

//                 float denom = Ix*Ix + Iy*Iy + 1e-4f;

//                 flow_x[i] = -Ix * It / denom;
//                 flow_y[i] = -Iy * It / denom;
//             }
//         }

//         //-------------------------------
//         // 3. Alle Punkte bewegen + zeichnen
//         //-------------------------------
//         for(int p = 0; p < pcount; p++)
//         {
//             int ix = (int)px[p];
//             int iy = (int)py[p];

//             if(ix < 1 || ix >= WIDTH-1 || iy < 1 || iy >= HEIGHT-1)
//                 continue;

//             int i = idx(ix, iy);

//             float vx = flow_x[i];
//             float vy = flow_y[i];

//             float nx = px[p] + vx * scale;
//             float ny = py[p] + vy * scale;

//             draw_line(flow_vis, (int)px[p], (int)py[p], (int)nx, (int)ny);

//             px[p] = nx;
//             py[p] = ny;
//         }
//     }

//     //------------------------------------------------------------------
//     // Ausgabe
//     //------------------------------------------------------------------
//     for(int i = 0; i < NPIX*3; i++)
//     {
//         printf("%u\n", flow_vis[i]);
//     }

//     return 0;
// }




//Calculation for the whole video just like before but drawing less lines as output


// #include <stdio.h>
// #include <stdint.h>
// #include <math.h>
// #include <stdlib.h>

// #include "../../../input/videos/split/split_video/headers/images_index.h"

// #define WIDTH IMG_WIDTH
// #define HEIGHT IMG_HEIGHT
// #define NPIX (WIDTH*HEIGHT)

// // Flow
// static float flow_x[NPIX];
// static float flow_y[NPIX];

// // Visualisierung
// static uint8_t flow_vis[NPIX*3];

// static inline int idx(int x,int y)
// {
//     return y*WIDTH + x;
// }

// //--------------------------------------------------
// void draw_line(uint8_t* img, int x0, int y0, int x1, int y1)
// {
//     int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
//     int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
//     int err = dx + dy;

//     while(1)
//     {
//         if(x0 >= 0 && x0 < WIDTH && y0 >= 0 && y0 < HEIGHT)
//         {
//             int i = (y0 * WIDTH + x0) * 3;
//             img[i+0] = 0;
//             img[i+1] = 255;
//             img[i+2] = 255;
//         }

//         if(x0 == x1 && y0 == y1) break;

//         int e2 = 2 * err;
//         if(e2 >= dy) { err += dy; x0 += sx; }
//         if(e2 <= dx) { err += dx; y0 += sy; }
//     }
// }

// int main()
// {
//     const uint8_t* img1;
//     const uint8_t* img2;

//     int step = 15;
//     float scale = 3.0f;

//     //------------------------------------------------------------------
//     // Hintergrund = letztes Frame
//     //------------------------------------------------------------------
//     const uint8_t* last = frames[IMG_COUNT-1];

//     for(int i = 0; i < NPIX; i++)
//     {
//         flow_vis[i*3+0] = last[i];
//         flow_vis[i*3+1] = last[i];
//         flow_vis[i*3+2] = last[i];
//     }

//     //------------------------------------------------------------------
//     // Über alle Frames
//     //------------------------------------------------------------------
//     for(int f = 0; f < IMG_COUNT - 1; f++)
//     {
//         img1 = frames[f];
//         img2 = frames[f+1];

//         //-------------------------------
//         // 1. Flow berechnen (komplett)
//         //-------------------------------
//         for(int y=1;y<HEIGHT-1;y++)
//         {
//             for(int x=1;x<WIDTH-1;x++)
//             {
//                 int i = idx(x,y);

//                 float Ix = img1[idx(x+1,y)] - img1[idx(x-1,y)];
//                 float Iy = img1[idx(x,y+1)] - img1[idx(x,y-1)];
//                 float It = img2[i] - img1[i];

//                 float denom = Ix*Ix + Iy*Iy + 1e-4f;

//                 flow_x[i] = -Ix * It / denom;
//                 flow_y[i] = -Iy * It / denom;
//             }
//         }

//         //-------------------------------
//         // 2. Für jedes Grid: EINE Linie
//         //-------------------------------
//         for(int y = 0; y < HEIGHT; y += step)
//         {
//             for(int x = 0; x < WIDTH; x += step)
//             {
//                 int i = idx(x,y);

//                 float vx = flow_x[i];
//                 float vy = flow_y[i];

//                 float mag = sqrtf(vx*vx + vy*vy);

//                 // kleine Bewegungen ignorieren
//                 if(mag < 0.3f) continue;

//                 int nx = (int)(x + vx * scale);
//                 int ny = (int)(y + vy * scale);

//                 draw_line(flow_vis, x, y, nx, ny);
//             }
//         }
//     }

//     //------------------------------------------------------------------
//     // Ausgabe
//     //------------------------------------------------------------------
//     for(int i = 0; i < NPIX*3; i++)
//     {
//         printf("%u\n", flow_vis[i]);
//     }

//     return 0;
// }




//Calculation for the whole video, calculate the average movement the video and draw a line for this
//Result can be seen in flow_lines_7.png

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

// Visualisierung
static uint8_t flow_vis[NPIX*3];

static inline int idx(int x,int y)
{
    return y*WIDTH + x;
}

//--------------------------------------------------
void draw_line(uint8_t* img, int x0, int y0, int x1, int y1)
{
    int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = dx + dy;

    while(1)
    {
        if(x0 >= 0 && x0 < WIDTH && y0 >= 0 && y0 < HEIGHT)
        {
            int i = (y0 * WIDTH + x0) * 3;
            img[i+0] = 255;   // rot für globalen Vektor
            img[i+1] = 0;
            img[i+2] = 0;
        }

        if(x0 == x1 && y0 == y1) break;

        int e2 = 2 * err;
        if(e2 >= dy) { err += dy; x0 += sx; }
        if(e2 <= dx) { err += dx; y0 += sy; }
    }
}

int main()
{
    const uint8_t* img1;
    const uint8_t* img2;

    //--------------------------------------------------
    // Akkumulatoren für globalen Flow
    //--------------------------------------------------
    double sum_vx = 0.0;
    double sum_vy = 0.0;
    long count = 0;

    float threshold = 0.3f;

    //--------------------------------------------------
    // Über alle Frames
    //--------------------------------------------------
    for(int f = 0; f < IMG_COUNT - 1; f++)
    {
        img1 = frames[f];
        img2 = frames[f+1];

        //-------------------------------
        // Flow berechnen
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
        // Mittelwert sammeln
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
    // Durchschnitt berechnen
    //--------------------------------------------------
    float avg_vx = 0.0f;
    float avg_vy = 0.0f;

    if(count > 0)
    {
        avg_vx = sum_vx / count;
        avg_vy = sum_vy / count;
    }

    printf("Average Flow: vx=%f vy=%f\n", avg_vx, avg_vy);

    //--------------------------------------------------
    // Hintergrund = letztes Frame
    //--------------------------------------------------
    const uint8_t* last = frames[IMG_COUNT-1];

    for(int i = 0; i < NPIX; i++)
    {
        flow_vis[i*3+0] = last[i];
        flow_vis[i*3+1] = last[i];
        flow_vis[i*3+2] = last[i];
    }

    //--------------------------------------------------
    // Globalen Vektor zeichnen (Bildmitte)
    //--------------------------------------------------
    int cx = WIDTH / 2;
    int cy = HEIGHT / 2;

    float scale = 50.0f; // größer für Sichtbarkeit

    int nx = (int)(cx + avg_vx * scale);
    int ny = (int)(cy + avg_vy * scale);

    draw_line(flow_vis, cx, cy, nx, ny);

    //--------------------------------------------------
    // Ausgabe
    //--------------------------------------------------
    for(int i = 0; i < NPIX*3; i++)
    {
        printf("%u\n", flow_vis[i]);
    }

    return 0;
}


