#include "./header/video_14/image_0000.h"
#include "./header/video_14/image_0001.h"
#define WIDTH 1280
#define HEIGHT 720
#include <math.h>
#include <stdlib.h>

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
    int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
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
    const uint8_t* img1 = image_0000;
    const uint8_t* img2 = image_0001;

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

            flow_x[i] = -Ix * It / denom;
            flow_y[i] = -Iy * It / denom;
        }
    }

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
    // 5. ARRAY AUSGABE (RGB)
    //------------------------------------------------------------------

    // for(int i = 0; i < NPIX*3; i++)
    // {
    //     printf("%u\n", flow_vis[i]);
    // }

    return 0;
}
