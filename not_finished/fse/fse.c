/******************************************************************************
 *
 *  Frequency Selective Extrapolation
 *  Teil 1
 *
 ******************************************************************************/

#include <stdint.h>
#include <stdbool.h>

#include "images/mask_image.h"

/******************************************************************************/
// Konfiguration
/******************************************************************************/

#define WIDTH 600
#define HEIGHT 600

#define MEDIAN_SIZE 9

#define ERROR_THRESHOLD 35

#define MAX_ITERATIONS 8

/******************************************************************************/
// Speicher
/******************************************************************************/

static uint8_t workImage[HEIGHT][WIDTH];

static uint8_t errorMask[HEIGHT][WIDTH];

static int16_t gradientX[HEIGHT][WIDTH];

static int16_t gradientY[HEIGHT][WIDTH];

/******************************************************************************/

static inline uint8_t absDiff(uint8_t a,uint8_t b)
{
    if(a>b)
        return a-b;

    return b-a;
}

/******************************************************************************/

static void copyImage(void)
{
    uint32_t x;
    uint32_t y;

    for(y=0;y<HEIGHT;y++)
    {
        for(x=0;x<WIDTH;x++)
        {
            workImage[y][x]=mask_image[y][x];
        }
    }
}

/******************************************************************************/

static void swap(uint8_t *a,uint8_t *b)
{
    uint8_t t=*a;
    *a=*b;
    *b=t;
}

/******************************************************************************/

static uint8_t median9(uint8_t v[9])
{
    uint8_t i;
    uint8_t j;

    for(i=0;i<8;i++)
    {
        for(j=i+1;j<9;j++)
        {
            if(v[j]<v[i])
            {
                swap(&v[i],&v[j]);
            }
        }
    }

    return v[4];
}

/******************************************************************************/

static uint8_t localMedian(uint32_t x,uint32_t y)
{
    uint8_t values[9];

    uint8_t k=0;

    int32_t xx;
    int32_t yy;

    for(yy=-1;yy<=1;yy++)
    {
        for(xx=-1;xx<=1;xx++)
        {
            values[k++]=workImage[y+yy][x+xx];
        }
    }

    return median9(values);
}

/******************************************************************************/

static void detectErrors(void)
{
    uint32_t x;
    uint32_t y;

    uint8_t med;

    for(y=0;y<HEIGHT;y++)
    {
        for(x=0;x<WIDTH;x++)
        {
            errorMask[y][x]=0;
        }
    }

    for(y=1;y<HEIGHT-1;y++)
    {
        for(x=1;x<WIDTH-1;x++)
        {
            med=localMedian(x,y);

            if(absDiff(workImage[y][x],med)>ERROR_THRESHOLD)
            {
                errorMask[y][x]=1;
            }
        }
    }
}

/******************************************************************************/

static uint32_t countErrors(void)
{
    uint32_t x;
    uint32_t y;

    uint32_t c=0;

    for(y=0;y<HEIGHT;y++)
    {
        for(x=0;x<WIDTH;x++)
        {
            c+=errorMask[y][x];
        }
    }

    return c;
}


/******************************************************************************
 *
 *  Teil 2
 *  Gradienten- und Richtungsberechnung
 *
 ******************************************************************************/

static inline int16_t abs16(int16_t v)
{
    return (v < 0) ? -v : v;
}

/******************************************************************************/

static void computeGradients(void)
{
    uint32_t x;
    uint32_t y;

    int16_t gx;
    int16_t gy;

    for(y = 0; y < HEIGHT; y++)
    {
        for(x = 0; x < WIDTH; x++)
        {
            gradientX[y][x] = 0;
            gradientY[y][x] = 0;
        }
    }

    for(y = 1; y < HEIGHT - 1; y++)
    {
        for(x = 1; x < WIDTH - 1; x++)
        {
            gx =
               -workImage[y-1][x-1]
               +workImage[y-1][x+1]
             -2*workImage[y][x-1]
             +2*workImage[y][x+1]
               -workImage[y+1][x-1]
               +workImage[y+1][x+1];

            gy =
               -workImage[y-1][x-1]
             -2*workImage[y-1][x]
               -workImage[y-1][x+1]
               +workImage[y+1][x-1]
             +2*workImage[y+1][x]
               +workImage[y+1][x+1];

            gradientX[y][x] = gx;
            gradientY[y][x] = gy;
        }
    }
}

/******************************************************************************/

typedef enum
{
    DIR_HORIZONTAL = 0,
    DIR_VERTICAL,
    DIR_DIAG1,
    DIR_DIAG2

} Direction;

/******************************************************************************/

static Direction dominantDirection(uint32_t x,uint32_t y)
{
    int16_t gx = gradientX[y][x];
    int16_t gy = gradientY[y][x];

    int16_t ax = abs16(gx);
    int16_t ay = abs16(gy);

    if(ax > (ay << 1))
    {
        return DIR_VERTICAL;
    }

    if(ay > (ax << 1))
    {
        return DIR_HORIZONTAL;
    }

    if((gx ^ gy) >= 0)
    {
        return DIR_DIAG1;
    }

    return DIR_DIAG2;
}

/******************************************************************************/

static uint8_t average2(uint8_t a,uint8_t b)
{
    return (uint8_t)(((uint16_t)a + (uint16_t)b) >> 1);
}

/******************************************************************************/

static uint8_t average4(uint8_t a,
                        uint8_t b,
                        uint8_t c,
                        uint8_t d)
{
    return (uint8_t)
    (
        (
            (uint16_t)a +
            (uint16_t)b +
            (uint16_t)c +
            (uint16_t)d
        ) >> 2
    );
}

/******************************************************************************/

static uint8_t extrapolateDirectional(uint32_t x,uint32_t y)
{
    switch(dominantDirection(x,y))
    {

        case DIR_HORIZONTAL:

            return average2(
                workImage[y][x-1],
                workImage[y][x+1]);

        case DIR_VERTICAL:

            return average2(
                workImage[y-1][x],
                workImage[y+1][x]);

        case DIR_DIAG1:

            return average2(
                workImage[y-1][x-1],
                workImage[y+1][x+1]);

        default:

            return average2(
                workImage[y-1][x+1],
                workImage[y+1][x-1]);
    }
}

/******************************************************************************
 *
 * Teil 3
 * Frequency Selective Extrapolation (iterativ)
 *
 ******************************************************************************/

static inline uint8_t isValidPixel(int32_t x,int32_t y)
{
    if(x < 0) return 0;
    if(y < 0) return 0;
    if(x >= WIDTH) return 0;
    if(y >= HEIGHT) return 0;

    return (errorMask[y][x] == 0);
}

/******************************************************************************/

static uint8_t directionalEstimate(uint32_t x,
                                   uint32_t y,
                                   Direction dir,
                                   uint8_t *valid)
{
    uint16_t sum = 0;
    uint8_t count = 0;

    *valid = 0;

    switch(dir)
    {
        case DIR_HORIZONTAL:

            if(isValidPixel(x-1,y))
            {
                sum += workImage[y][x-1];
                count++;
            }

            if(isValidPixel(x+1,y))
            {
                sum += workImage[y][x+1];
                count++;
            }

            break;

        case DIR_VERTICAL:

            if(isValidPixel(x,y-1))
            {
                sum += workImage[y-1][x];
                count++;
            }

            if(isValidPixel(x,y+1))
            {
                sum += workImage[y+1][x];
                count++;
            }

            break;

        case DIR_DIAG1:

            if(isValidPixel(x-1,y-1))
            {
                sum += workImage[y-1][x-1];
                count++;
            }

            if(isValidPixel(x+1,y+1))
            {
                sum += workImage[y+1][x+1];
                count++;
            }

            break;

        default:

            if(isValidPixel(x-1,y+1))
            {
                sum += workImage[y+1][x-1];
                count++;
            }

            if(isValidPixel(x+1,y-1))
            {
                sum += workImage[y-1][x+1];
                count++;
            }

            break;
    }

    if(count == 0)
        return 0;

    *valid = 1;

    return (uint8_t)(sum / count);
}

/******************************************************************************/

static uint8_t reconstructPixel(uint32_t x,uint32_t y)
{
    uint8_t estimate[4];
    uint8_t valid[4];

    uint16_t sum = 0;
    uint8_t count = 0;
    uint8_t i;

    for(i=0;i<4;i++)
    {
        estimate[i] =
            directionalEstimate(x,y,(Direction)i,&valid[i]);

        if(valid[i])
        {
            sum += estimate[i];
            count++;
        }
    }

    if(count == 0)
    {
        return localMedian(x,y);
    }

    return (uint8_t)(sum / count);
}

/******************************************************************************/

static void reconstructIteration(void)
{
    uint32_t x;
    uint32_t y;

    uint8_t changed = 0;

    for(y=1;y<HEIGHT-1;y++)
    {
        for(x=1;x<WIDTH-1;x++)
        {
            if(errorMask[y][x])
            {
                workImage[y][x] = reconstructPixel(x,y);

                errorMask[y][x] = 0;

                changed = 1;
            }
        }
    }

    (void)changed;
}

/******************************************************************************/

static void reconstructImage(void)
{
    uint32_t i;

    for(i=0;i<MAX_ITERATIONS;i++)
    {
        computeGradients();

        reconstructIteration();
    }
}

/******************************************************************************
 *
 * Teil 4
 * Richtungsgewichtung
 *
 ******************************************************************************/

typedef struct
{
    uint8_t value;
    uint16_t energy;
    uint8_t valid;

} DirectionEstimate;

/******************************************************************************/

static uint16_t directionEnergy(uint32_t x,
                                uint32_t y,
                                Direction dir)
{
    int16_t g1;
    int16_t g2;

    switch(dir)
    {
        case DIR_HORIZONTAL:

            g1 = gradientX[y][x-1];
            g2 = gradientX[y][x+1];
            break;

        case DIR_VERTICAL:

            g1 = gradientY[y-1][x];
            g2 = gradientY[y+1][x];
            break;

        case DIR_DIAG1:

            g1 = gradientX[y-1][x-1] +
                 gradientY[y-1][x-1];

            g2 = gradientX[y+1][x+1] +
                 gradientY[y+1][x+1];
            break;

        default:

            g1 = gradientX[y+1][x-1] -
                 gradientY[y+1][x-1];

            g2 = gradientX[y-1][x+1] -
                 gradientY[y-1][x+1];
            break;
    }

    return (uint16_t)(abs16(g1) + abs16(g2));
}

/******************************************************************************/

static void estimateDirection(uint32_t x,
                              uint32_t y,
                              Direction dir,
                              DirectionEstimate *est)
{
    est->value =
        directionalEstimate(x,y,dir,&est->valid);

    if(est->valid)
    {
        est->energy = directionEnergy(x,y,dir) + 1;
    }
    else
    {
        est->energy = 0;
    }
}

/******************************************************************************/

static uint8_t weightedEstimate(uint32_t x,uint32_t y)
{
    DirectionEstimate est[4];

    uint32_t weightedSum = 0;
    uint32_t totalWeight = 0;

    uint8_t i;

    for(i=0;i<4;i++)
    {
        estimateDirection(x,y,(Direction)i,&est[i]);

        if(est[i].valid)
        {
            weightedSum +=
                ((uint32_t)est[i].value) *
                ((uint32_t)est[i].energy);

            totalWeight += est[i].energy;
        }
    }

    if(totalWeight == 0)
        return localMedian(x,y);

    return (uint8_t)(weightedSum / totalWeight);
}

/******************************************************************************
 *
 * Teil 5
 * Adaptive Suche entlang der Frequenzrichtung
 *
 ******************************************************************************/

#define SEARCH_RADIUS 7

/*****************************************************************************/

static uint8_t searchDirection(uint32_t x,
                               uint32_t y,
                               int dx,
                               int dy,
                               uint8_t *found)
{
    int32_t xx;
    int32_t yy;

    uint32_t r;

    *found = 0;

    for(r = 1; r <= SEARCH_RADIUS; r++)
    {
        xx = (int32_t)x + dx * (int32_t)r;
        yy = (int32_t)y + dy * (int32_t)r;

        if(xx < 0)
            break;

        if(yy < 0)
            break;

        if(xx >= WIDTH)
            break;

        if(yy >= HEIGHT)
            break;

        if(errorMask[yy][xx] == 0)
        {
            *found = 1;
            return workImage[yy][xx];
        }
    }

    return 0;
}


/*****************************************************************************/

static uint8_t extrapolateAdaptive(uint32_t x,uint32_t y)
{
    uint8_t left;
    uint8_t right;

    uint8_t ok1;
    uint8_t ok2;

    Direction dir;

    dir = dominantDirection(x,y);

    switch(dir)
    {

        case DIR_HORIZONTAL:

            left  = searchDirection(x,y,-1,0,&ok1);
            right = searchDirection(x,y, 1,0,&ok2);

            break;

        case DIR_VERTICAL:

            left  = searchDirection(x,y,0,-1,&ok1);
            right = searchDirection(x,y,0, 1,&ok2);

            break;

        case DIR_DIAG1:

            left  = searchDirection(x,y,-1,-1,&ok1);
            right = searchDirection(x,y, 1, 1,&ok2);

            break;

        default:

            left  = searchDirection(x,y,-1, 1,&ok1);
            right = searchDirection(x,y, 1,-1,&ok2);

            break;
    }

    if(ok1 && ok2)
    {
        return (uint8_t)(((uint16_t)left + right) >> 1);
    }

    if(ok1)
        return left;

    if(ok2)
        return right;

    return localMedian(x,y);
}

int main(void)
{
    copyImage();

    detectErrors();

    reconstructImage();

    for(uint32_t y=0;y<HEIGHT;y++){
    for(uint32_t x=0;x<WIDTH;x++)
        {
            printf("%3u\n",workImage[y][x]);
        }

        //printf("\n");
}
   
    return 0;
}