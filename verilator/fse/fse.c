

#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include "images/image_data.h"

#define MAX_FREQS 32
#define LUT_SIZE 1024

#define PI 3.14159265358979323846f

#define IMAGE_WIDTH 600
#define IMAGE_HEIGHT 600

#define MAX_IMAGE_WIDTH 600
#define MAX_IMAGE_HEIGHT 600

/*
    Marker für unbekannte Pixel.
    Falls deine Bilddaten einen anderen Wert verwenden:
*/
#define UNKNOWN_PIXEL 0.0f

typedef struct
{
    int16_t u;
    int16_t v;

    float amplitude;
    float phase;

} FSE_Component;

static float residual[MAX_IMAGE_HEIGHT][MAX_IMAGE_WIDTH];

static int output_image[MAX_IMAGE_HEIGHT][MAX_IMAGE_WIDTH];

static FSE_Component components[MAX_FREQS];

static float sin_lut[LUT_SIZE];
static float cos_lut[LUT_SIZE];

static inline int lut_index(float angle)
{
    while(angle < 0.0f)
        angle += 2.0f * PI;

    while(angle >= 2.0f * PI)
        angle -= 2.0f * PI;

    int idx =
        (int)(angle * LUT_SIZE / (2.0f * PI));

    if(idx >= LUT_SIZE)
        idx = LUT_SIZE - 1;

    return idx;
}

static inline float fast_cos(float angle)
{
    return cos_lut[lut_index(angle)];
}

void fse_init(void)
{
    for(int i = 0; i < LUT_SIZE; i++)
    {
        float a =
            2.0f * PI *
            (float)i /
            (float)LUT_SIZE;

        sin_lut[i] = sinf(a);
        cos_lut[i] = cosf(a);
    }
}

static void copy_source_to_residual(void)
{
    for(int y = 0; y < IMAGE_HEIGHT; y++)
    {
        for(int x = 0; x < IMAGE_WIDTH; x++)
        {
            residual[y][x] =
                image_data[y][x];
        }
    }
}

static void copy_to_output(void)
{
    for(int y = 0; y < IMAGE_HEIGHT; y++)
    {
        for(int x = 0; x < IMAGE_WIDTH; x++)
        {
            output_image[y][x] =
                (int)image_data[y][x];
        }
    }
}

/*
    Sehr primitive Frequenzschätzung
    (nur zwei Schleifen)
*/
static void find_dominant_frequency(
    float img[MAX_IMAGE_HEIGHT][MAX_IMAGE_WIDTH],
    FSE_Component* comp)
{
    float dx = 0.0f;
    float dy = 0.0f;

    uint32_t count = 0;

    for(int y = 0; y < IMAGE_HEIGHT - 1; y++)
    {
        for(int x = 0; x < IMAGE_WIDTH - 1; x++)
        {
            if(img[y][x] == UNKNOWN_PIXEL)
                continue;

            if(img[y][x + 1] != UNKNOWN_PIXEL)
            {
                dx += fabsf(
                    img[y][x + 1] -
                    img[y][x]);
            }

            if(img[y + 1][x] != UNKNOWN_PIXEL)
            {
                dy += fabsf(
                    img[y + 1][x] -
                    img[y][x]);
            }

            count++;
        }
    }

    if(dx > dy)
    {
        comp->u = IMAGE_WIDTH / 8;
        comp->v = 0;
    }
    else
    {
        comp->u = 0;
        comp->v = IMAGE_HEIGHT / 8;
    }

    comp->phase = 0.0f;

    if(count == 0)
        count = 1;

    comp->amplitude =
        (dx + dy) /
        (float)count;
}

static void subtract_component(
    float img[MAX_IMAGE_HEIGHT][MAX_IMAGE_WIDTH],
    const FSE_Component* c)
{
    for(int y = 0; y < IMAGE_HEIGHT; y++)
    {
        for(int x = 0; x < IMAGE_WIDTH; x++)
        {
            if(img[y][x] == UNKNOWN_PIXEL)
                continue;

            float angle =
                2.0f * PI *
                (
                    ((float)c->u * x /
                     IMAGE_WIDTH)
                    +
                    ((float)c->v * y /
                     IMAGE_HEIGHT)
                );

            img[y][x] -=
                c->amplitude *
                fast_cos(angle);
        }
    }
}

void fse_analyse(int component_count)
{
    if(component_count > MAX_FREQS)
        component_count = MAX_FREQS;

    copy_source_to_residual();

    for(int i = 0; i < component_count; i++)
    {
        find_dominant_frequency(
            residual,
            &components[i]);

        subtract_component(
            residual,
            &components[i]);
    }
}

void fill_unknown_pixels(int component_count)
{
    copy_to_output();

    if(component_count > MAX_FREQS)
        component_count = MAX_FREQS;

    for(int y = 0; y < IMAGE_HEIGHT; y++)
    {
        for(int x = 0; x < IMAGE_WIDTH; x++)
        {
            if(output_image[y][x] != UNKNOWN_PIXEL)
                continue;

            float value = 0.0f;

            for(int k = 0; k < component_count; k++)
            {
                float angle =
                    2.0f * PI *
                    (
                        ((float)components[k].u * x /
                         IMAGE_WIDTH)
                        +
                        ((float)components[k].v * y /
                         IMAGE_HEIGHT)
                    );

                value +=
                    components[k].amplitude *
                    fast_cos(angle);
            }

            output_image[y][x] = (int)value;
        }
    }
}

int main(void)
{
    fse_init();

    fse_analyse(16);

    fill_unknown_pixels(16);

    for(int y = 0; y < IMAGE_HEIGHT; y++)
    {
        for(int x = 0; x < IMAGE_WIDTH; x++)
        {
            printf("%i\n", output_image[y][x]);
        }
    }

    return 0;
}