


//The following code draws lines in the output image depicting the actual movement

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include "uart_com.h"


#include "../input/videos/split/split_video/headers/image_0000.h"
#include "../input/videos/split/split_video/headers/image_0001.h"

#define WIDTH 426
#define HEIGHT 240
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
// Bresenham Linie
//--------------------------------------------------
void draw_line(uint8_t* img, int x0, int y0, int x1, int y1)
{
    int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1; //!MATH
    int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1; //!MATH
    int err = dx + dy;

    while(1)
    {
        if(x0 >= 0 && x0 < WIDTH && y0 >= 0 && y0 < HEIGHT)
        {
            int i = (y0 * WIDTH + x0) * 3;

            // Farbe (cyan)
            img[i+0] = 0;
            img[i+1] = 255;
            img[i+2] = 255;
        }

        if(x0 == x1 && y0 == y1) break;

        int e2 = 2 * err;
        if(e2 >= dy) { err += dy; x0 += sx; }
        if(e2 <= dx) { err += dx; y0 += sy; }
    }
}

int main()
{

    uart_puts("Start\n");
    const uint8_t* img1 = image_0000;
    const uint8_t* img2 = image_0001;

    //------------------------------------------------------------------
    // 1. Optical Flow berechnen
    //------------------------------------------------------------------

    for(int y=1;y<HEIGHT-1;y++)
    {
        for(int x=1;x<WIDTH-1;x++)
        {
            uart_puts("In Loop\n");
            int i = idx(x,y);
            uart_puts("In Loop 2\n");
            // float Ix = img1[idx(x+1,y)] - img1[idx(x-1,y)];
            // float Iy = img1[idx(x,y+1)] - img1[idx(x,y-1)];
            // float It = img2[i] - img1[i];
            float Ix = 4.0f;
            float Iy = 2.0f;
            float It = 6.0f;
            uart_puts("In Loop\n");
            float denom = Ix*Ix + Iy*Iy + 1e-4f;

            flow_x[i] = -Ix * It / denom;
            flow_y[i] = -Iy * It / denom;
        }
    }

    uart_puts("Oprflow Calc ready");

    //------------------------------------------------------------------
    // 2. Flow glätten (3x3 Filter)
    //------------------------------------------------------------------

    for(int y=1;y<HEIGHT-1;y++)
    {
        for(int x=1;x<WIDTH-1;x++)
        {
            float sumx = 0, sumy = 0;

            for(int dy=-1; dy<=1; dy++)
            for(int dx=-1; dx<=1; dx++)
            {
                int j = idx(x+dx, y+dy);
                sumx += flow_x[j];
                sumy += flow_y[j];
            }

            int i = idx(x,y);
            flow_x[i] = sumx / 9.0f;
            flow_y[i] = sumy / 9.0f;
        }
    }

    uart_puts("Next Step done");
    //------------------------------------------------------------------
    // 3. Hintergrund setzen (Graubild)
    //------------------------------------------------------------------

    for(int i = 0; i < NPIX; i++)
    {
        flow_vis[i*3+0] = img1[i];
        flow_vis[i*3+1] = img1[i];
        flow_vis[i*3+2] = img1[i];
    }

    //------------------------------------------------------------------
    // 4. Trajektorien zeichnen
    //------------------------------------------------------------------

    int step = 15;      // dichter als vorher
    float scale = 3.0f;
    int steps = 10;     // Länge der Linien

    for(int y = 0; y < HEIGHT; y += step)
    {
        for(int x = 0; x < WIDTH; x += step)
        {
            uart_puts("Drawing");
            int i = idx(x,y);

            float vx0 = flow_x[i];
            float vy0 = flow_y[i];

            float mag = sqrtf(vx0*vx0 + vy0*vy0); //!MATH

            // nur signifikante Bewegung
            if(mag < 0.3f) continue;

            float px = x;
            float py = y;

            for(int s = 0; s < steps; s++)
            {
                int ix = (int)px;
                int iy = (int)py;

                if(ix < 0 || ix >= WIDTH || iy < 0 || iy >= HEIGHT)
                    break;

                int j = idx(ix, iy);

                float vx = flow_x[j];
                float vy = flow_y[j];

                float nx = px + vx * scale;
                float ny = py + vy * scale;

                draw_line(flow_vis, (int)px, (int)py, (int)nx, (int)ny);

                px = nx;
                py = ny;
            }
        }
    }

    uart_puts("Done");
    //------------------------------------------------------------------
    // 5. ARRAY AUSGABE (RGB)
    //------------------------------------------------------------------

    for(int i = 0; i < NPIX*3; i++)
    {
        uart_putint(flow_vis[i]);
        uart_puts("\n");
    }

    return 0;
}
