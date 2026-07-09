//Result is result4.pgm

#include <stdio.h>
#include <stdint.h>
#include <math.h>


//#include "fse.h"
#include "images/mask_image.h"

#define IMAGE_WIDTH 600
#define IMAGE_HEIGHT 600

#define MAX_IMAGE_WIDTH 600
#define MAX_IMAGE_HEIGHT 600

#define MAX_FREQS 32

#define LUT_SIZE 2048

#define PI 3.14159265358979323846f

typedef struct
{
    int16_t u;
    int16_t v;

    float amplitude;
    float phase;

} FSE_Component;

static float work_image[MAX_IMAGE_HEIGHT][MAX_IMAGE_WIDTH];

static uint8_t unknown_mask[MAX_IMAGE_HEIGHT][MAX_IMAGE_WIDTH];

static uint8_t output_image[MAX_IMAGE_HEIGHT][MAX_IMAGE_WIDTH];

static FSE_Component components[MAX_FREQS];

static float sin_lut[LUT_SIZE];
static float cos_lut[LUT_SIZE];



static inline int lut_index(float angle)
{
    while(angle < 0.0f)
        angle += 2.0f * PI;

    while(angle >= 2.0f * PI)
        angle -= 2.0f * PI;

    int index =
        (int)(angle * LUT_SIZE / (2.0f * PI));

    if(index >= LUT_SIZE)
        index = LUT_SIZE - 1;

    return index;
}

static inline float fast_sin(float angle)
{
    return sin_lut[lut_index(angle)];
}

static inline float fast_cos(float angle)
{
    return cos_lut[lut_index(angle)];
}

void fse_init(void)
{
    for(int i=0;i<LUT_SIZE;i++)
    {
        float a =
            2.0f *
            PI *
            (float)i /
            (float)LUT_SIZE;

        sin_lut[i] = sinf(a);
        cos_lut[i] = cosf(a);
    }
}

static void prepare_images(void)
{
    for(int y=0;y<IMAGE_HEIGHT;y++)
    {
        for(int x=0;x<IMAGE_WIDTH;x++)
        {
            uint8_t value =
                mask_image[y][x];

            output_image[y][x] = value;

            work_image[y][x] =
                (float)value;

            if(value <= 1)
            {
                unknown_mask[y][x] = 1;
            }
            else
            {
                unknown_mask[y][x] = 0;
            }
        }
    }
}

static inline int valid_pixel(int x,int y)
{
    return unknown_mask[y][x] == 0;
}

static inline uint8_t clamp(float value)
{
    if(value < 0.0f)
        value = 0.0f;

    if(value > 255.0f)
        value = 255.0f;

    return (uint8_t)(value + 0.5f);
}


void fse_print_image(void)
{
    for(int y=0;y<IMAGE_HEIGHT;y++)
    {
        for(int x=0;x<IMAGE_WIDTH;x++)
        {
            printf("%u\n",
                   output_image[y][x]);
        }
    }
}


/* ---------------------------------------------------------
   Fourier-Korrelation einer einzelnen Frequenz
   ---------------------------------------------------------*/

static void evaluate_frequency(
    int u,
    int v,
    FSE_Component *component)
{
    float re = 0.0f;
    float im = 0.0f;

    int samples = 0;

    for(int y=0;y<IMAGE_HEIGHT;y++)
    {
        for(int x=0;x<IMAGE_WIDTH;x++)
        {
            if(!valid_pixel(x,y))
                continue;

            float value = work_image[y][x];

            float angle =
                2.0f * PI *
                (
                    ((float)u * x / IMAGE_WIDTH) +
                    ((float)v * y / IMAGE_HEIGHT)
                );

            re += value * fast_cos(angle);
            im -= value * fast_sin(angle);

            samples++;
        }
    }

    if(samples == 0)
        samples = 1;

    component->u = u;
    component->v = v;

    component->amplitude =
        (2.0f * sqrtf(re * re + im * im))
        / (float)samples;

    component->phase = atan2f(im,re);
}



/* ---------------------------------------------------------
   Dominante Frequenzen bestimmen
   ---------------------------------------------------------*/

static void analyse_frequencies(void)
{
    int idx = 0;

    for(int v=0; v<16 && idx<MAX_FREQS; v++)
    {
        for(int u=0; u<16 && idx<MAX_FREQS; u++)
        {
            evaluate_frequency(u, v, &components[idx]);

            idx++;
        }
    }

    /* --- Sortiere nach Energie (Amplitude) --- */
    for(int i=0;i<idx-1;i++)
    {
        for(int j=i+1;j<idx;j++)
        {
            if(components[j].amplitude >
               components[i].amplitude)
            {
                FSE_Component tmp = components[i];
                components[i] = components[j];
                components[j] = tmp;
            }
        }
    }
}

void fse_process(void)
{
    prepare_images();

    analyse_frequencies();

    /* -------------------------------------------------
       1. Frequenzen nach Stärke sortieren (WICHTIG!)
       ------------------------------------------------- */
    for(int i = 0; i < MAX_FREQS - 1; i++)
    {
        for(int j = i + 1; j < MAX_FREQS; j++)
        {
            if(components[j].amplitude > components[i].amplitude)
            {
                FSE_Component tmp = components[i];
                components[i] = components[j];
                components[j] = tmp;
            }
        }
    }

    /* -------------------------------------------------
       2. Rekonstruktion NUR beschädigter Pixel
       ------------------------------------------------- */
    for(int y = 0; y < IMAGE_HEIGHT; y++)
    {
        for(int x = 0; x < IMAGE_WIDTH; x++)
        {
            if(!unknown_mask[y][x])
            {
                output_image[y][x] = (uint8_t)work_image[y][x];
                continue;
            }

            float value = 0.0f;

            /* nur stärkste Frequenzen nutzen */
            for(int k = 0; k < MAX_FREQS / 2; k++)
            {
                float u = components[k].u;
                float v = components[k].v;

                float angle =
                    2.0f * PI *
                    (
                        (u * x / IMAGE_WIDTH) +
                        (v * y / IMAGE_HEIGHT)
                    );

                /* leichte Hochfrequenz-Dämpfung */
                float weight =
                    1.0f / (1.0f + 0.005f * (u*u + v*v));

                value +=
                    weight *
                    components[k].amplitude *
                    cosf(angle + components[k].phase);
            }

            /* konservativer Blend → weniger Veränderung */
            float original = work_image[y][x];

            float blended =
                0.60f * value +
                0.40f * original;

            output_image[y][x] =
                (uint8_t)(
                    blended < 0 ? 0 :
                    blended > 255 ? 255 :
                    blended + 0.5f
                );
        }
    }

    /* -------------------------------------------------
       3. sehr vorsichtige Stabilisierung
       (nur wenn wirklich nötig)
       ------------------------------------------------- */
    for(int iter = 0; iter < 2; iter++)
    {
        for(int y = 1; y < IMAGE_HEIGHT - 1; y++)
        {
            for(int x = 1; x < IMAGE_WIDTH - 1; x++)
            {
                if(!unknown_mask[y][x])
                    continue;

                float center = output_image[y][x];

                float avg =
                    center * 0.6f +
                    output_image[y - 1][x] * 0.1f +
                    output_image[y + 1][x] * 0.1f +
                    output_image[y][x - 1] * 0.1f +
                    output_image[y][x + 1] * 0.1f;

                output_image[y][x] = (uint8_t)(avg + 0.5f);
            }
        }
    }
}

int main(void)
{
    fse_init();

    fse_process();

    fse_print_image();

    return 0;
}

