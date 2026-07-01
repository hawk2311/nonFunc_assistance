

// #include <stdint.h>
// #include <stdio.h>
// #include <math.h>

// #include "images/image_data.h"

// #define MAX_FREQS 32
// #define LUT_SIZE 1024

// #define PI 3.14159265358979323846f

// #define IMAGE_WIDTH 600
// #define IMAGE_HEIGHT 600

// #define MAX_IMAGE_WIDTH 600
// #define MAX_IMAGE_HEIGHT 600

// /*
//     Marker für unbekannte Pixel.
//     Falls deine Bilddaten einen anderen Wert verwenden:
// */
// #define UNKNOWN_PIXEL 0.0f

// typedef struct
// {
//     int16_t u;
//     int16_t v;

//     float amplitude;
//     float phase;

// } FSE_Component;

// static float residual[MAX_IMAGE_HEIGHT][MAX_IMAGE_WIDTH];

// static int output_image[MAX_IMAGE_HEIGHT][MAX_IMAGE_WIDTH];

// static FSE_Component components[MAX_FREQS];

// static float sin_lut[LUT_SIZE];
// static float cos_lut[LUT_SIZE];

// static inline int lut_index(float angle)
// {
//     while(angle < 0.0f)
//         angle += 2.0f * PI;

//     while(angle >= 2.0f * PI)
//         angle -= 2.0f * PI;

//     int idx =
//         (int)(angle * LUT_SIZE / (2.0f * PI));

//     if(idx >= LUT_SIZE)
//         idx = LUT_SIZE - 1;

//     return idx;
// }

// static inline float fast_cos(float angle)
// {
//     return cos_lut[lut_index(angle)];
// }

// void fse_init(void)
// {
//     for(int i = 0; i < LUT_SIZE; i++)
//     {
//         float a =
//             2.0f * PI *
//             (float)i /
//             (float)LUT_SIZE;

//         sin_lut[i] = sinf(a);
//         cos_lut[i] = cosf(a);
//     }
// }

// static void copy_source_to_residual(void)
// {
//     for(int y = 0; y < IMAGE_HEIGHT; y++)
//     {
//         for(int x = 0; x < IMAGE_WIDTH; x++)
//         {
//             residual[y][x] =
//                 image_data[y][x];
//         }
//     }
// }

// static void copy_to_output(void)
// {
//     for(int y = 0; y < IMAGE_HEIGHT; y++)
//     {
//         for(int x = 0; x < IMAGE_WIDTH; x++)
//         {
//             output_image[y][x] =
//                 (int)image_data[y][x];
//         }
//     }
// }

// /*
//     Sehr primitive Frequenzschätzung
//     (nur zwei Schleifen)
// */
// static void find_dominant_frequency(
//     float img[MAX_IMAGE_HEIGHT][MAX_IMAGE_WIDTH],
//     FSE_Component* comp)
// {
//     float dx = 0.0f;
//     float dy = 0.0f;

//     uint32_t count = 0;

//     for(int y = 0; y < IMAGE_HEIGHT - 1; y++)
//     {
//         for(int x = 0; x < IMAGE_WIDTH - 1; x++)
//         {
//             if(img[y][x] == UNKNOWN_PIXEL)
//                 continue;

//             if(img[y][x + 1] != UNKNOWN_PIXEL)
//             {
//                 dx += fabsf(
//                     img[y][x + 1] -
//                     img[y][x]);
//             }

//             if(img[y + 1][x] != UNKNOWN_PIXEL)
//             {
//                 dy += fabsf(
//                     img[y + 1][x] -
//                     img[y][x]);
//             }

//             count++;
//         }
//     }

//     if(dx > dy)
//     {
//         comp->u = IMAGE_WIDTH / 8;
//         comp->v = 0;
//     }
//     else
//     {
//         comp->u = 0;
//         comp->v = IMAGE_HEIGHT / 8;
//     }

//     comp->phase = 0.0f;

//     if(count == 0)
//         count = 1;

//     comp->amplitude =
//         (dx + dy) /
//         (float)count;
// }

// static void subtract_component(
//     float img[MAX_IMAGE_HEIGHT][MAX_IMAGE_WIDTH],
//     const FSE_Component* c)
// {
//     for(int y = 0; y < IMAGE_HEIGHT; y++)
//     {
//         for(int x = 0; x < IMAGE_WIDTH; x++)
//         {
//             if(img[y][x] == UNKNOWN_PIXEL)
//                 continue;

//             float angle =
//                 2.0f * PI *
//                 (
//                     ((float)c->u * x /
//                      IMAGE_WIDTH)
//                     +
//                     ((float)c->v * y /
//                      IMAGE_HEIGHT)
//                 );

//             img[y][x] -=
//                 c->amplitude *
//                 fast_cos(angle);
//         }
//     }
// }

// void fse_analyse(int component_count)
// {
//     if(component_count > MAX_FREQS)
//         component_count = MAX_FREQS;

//     copy_source_to_residual();

//     for(int i = 0; i < component_count; i++)
//     {
//         find_dominant_frequency(
//             residual,
//             &components[i]);

//         subtract_component(
//             residual,
//             &components[i]);
//     }
// }

// void fill_unknown_pixels(int component_count)
// {
//     copy_to_output();

//     if(component_count > MAX_FREQS)
//         component_count = MAX_FREQS;

//     for(int y = 0; y < IMAGE_HEIGHT; y++)
//     {
//         for(int x = 0; x < IMAGE_WIDTH; x++)
//         {
//             if(output_image[y][x] != UNKNOWN_PIXEL)
//                 continue;

//             float value = 0.0f;

//             for(int k = 0; k < component_count; k++)
//             {
//                 float angle =
//                     2.0f * PI *
//                     (
//                         ((float)components[k].u * x /
//                          IMAGE_WIDTH)
//                         +
//                         ((float)components[k].v * y /
//                          IMAGE_HEIGHT)
//                     );

//                 value +=
//                     components[k].amplitude *
//                     fast_cos(angle);
//             }

//             output_image[y][x] = (int)value;
//         }
//     }
// }

// int main(void)
// {
//     fse_init();

//     fse_analyse(16);

//     fill_unknown_pixels(16);

//     for(int y = 0; y < IMAGE_HEIGHT; y++)
//     {
//         for(int x = 0; x < IMAGE_WIDTH; x++)
//         {
//             printf("%i\n", output_image[y][x]);
//         }
//     }

//     return 0;
// }


#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include "images/mask_image.h"

#define IMAGE_WIDTH      600
#define IMAGE_HEIGHT     600

#define MAX_FREQS        32

#define LUT_SIZE         1024

#define PI               3.14159265358979323846f

#define UNKNOWN_PIXEL    0.0f

typedef struct
{
    uint16_t u;
    uint16_t v;

    float amplitude;
    float phase;

}FSE_Component;


/*------------------------------------------------------------------
    Speicher
------------------------------------------------------------------*/

static float residual[IMAGE_HEIGHT][IMAGE_WIDTH];

static int output_image[IMAGE_HEIGHT][IMAGE_WIDTH];

static FSE_Component components[MAX_FREQS];

static float sin_lut[LUT_SIZE];

static float cos_lut[LUT_SIZE];


/*------------------------------------------------------------------
    LUT
------------------------------------------------------------------*/

static inline int lut_index(float angle)
{
    while(angle < 0.0f)
        angle += 2.0f*PI;

    while(angle >= 2.0f*PI)
        angle -= 2.0f*PI;

    int idx=(int)(angle*LUT_SIZE/(2.0f*PI));

    if(idx>=LUT_SIZE)
        idx=LUT_SIZE-1;

    return idx;
}

static inline float fast_cos(float angle)
{
    return cos_lut[lut_index(angle)];
}

static inline float fast_sin(float angle)
{
    return sin_lut[lut_index(angle)];
}

void fse_init(void)
{
    for(int i=0;i<LUT_SIZE;i++)
    {
        float a=
            2.0f*
            PI*
            i/
            LUT_SIZE;

        sin_lut[i]=sinf(a);

        cos_lut[i]=cosf(a);
    }
}


/*------------------------------------------------------------------
    Bild kopieren
------------------------------------------------------------------*/

// static void copy_input(void)
// {
//     for(int y=0;y<IMAGE_HEIGHT;y++)
//     {
//         for(int x=0;x<IMAGE_WIDTH;x++)
//         {
//             residual[y][x]=
//                 mask_image[y][x];

//             output_image[y][x]=
//                 mask_image[y][x];
//         }
//     }
// }


/*------------------------------------------------------------------
    Frequenzanalyse vorbereiten
------------------------------------------------------------------*/

static float mean_value(void)
{
    float sum=0.0f;

    uint32_t count=0;

    for(int y=0;y<IMAGE_HEIGHT;y++)
    {
        for(int x=0;x<IMAGE_WIDTH;x++)
        {
            if(residual[y][x]==UNKNOWN_PIXEL)
                continue;

            sum+=residual[y][x];

            count++;
        }
    }

    if(count==0)
        return 0.0f;

    return sum/(float)count;
}


static void find_dominant_frequency(
    float img[IMAGE_HEIGHT][IMAGE_WIDTH],
    FSE_Component* comp)
{
    float best_energy = 0.0f;

    int best_u = 0;
    int best_v = 0;

    /*
        Sehr einfache Projektion:
        Wir testen nur wenige "Richtungen",
        statt kompletten Frequenzraum.
    */

    for(int v = 0; v < IMAGE_HEIGHT/2; v += 8)
    {
        for(int u = 0; u < IMAGE_WIDTH/2; u += 8)
        {
            float sum = 0.0f;

            for(int y = 0; y < IMAGE_HEIGHT; y += 4)
            {
                for(int x = 0; x < IMAGE_WIDTH; x += 4)
                {
                    float p = img[y][x];

                    if(p == UNKNOWN_PIXEL)
                        continue;

                    float angle =
                        2.0f * PI *
                        (
                            ((float)u * x / IMAGE_WIDTH) +
                            ((float)v * y / IMAGE_HEIGHT)
                        );

                    sum += p * fast_cos(angle);
                }
            }

            float energy = sum * sum;

            if(energy > best_energy)
            {
                best_energy = energy;
                best_u = u;
                best_v = v;
            }
        }
    }

    comp->u = best_u;
    comp->v = best_v;
    comp->phase = 0.0f;
    comp->amplitude = sqrtf(best_energy) / (IMAGE_WIDTH * IMAGE_HEIGHT);
}


void fse_analyse(int component_count)
{
    if(component_count > MAX_FREQS)
        component_count = MAX_FREQS;

    //copy_input();

    float m = mean_value();

    /*
        Residual initialisieren:
        unbekannte Pixel werden auf Mittelwert gesetzt
        -> stabilisiert Frequenzsuche
    */
    for(int y=0;y<IMAGE_HEIGHT;y++)
    {
        for(int x=0;x<IMAGE_WIDTH;x++)
        {
            if(residual[y][x] == UNKNOWN_PIXEL)
                residual[y][x] = m;
        }
    }

    for(int i=0;i<component_count;i++)
    {
        find_dominant_frequency(residual, &components[i]);

        /*
            Sehr einfache Subtraktion (kein echtes Orthogonal Matching,
            aber stabil genug für Embedded)
        */

        for(int y=0;y<IMAGE_HEIGHT;y++)
        {
            for(int x=0;x<IMAGE_WIDTH;x++)
            {
                float angle =
                    2.0f * PI *
                    (
                        ((float)components[i].u * x / IMAGE_WIDTH) +
                        ((float)components[i].v * y / IMAGE_HEIGHT)
                    );

                residual[y][x] -=
                    components[i].amplitude *
                    fast_cos(angle);
            }
        }
    }
}


void fse_reconstruct_and_fill(void)
{
    /*
        Start: Originalbild kopieren
    */
    for(int y=0;y<IMAGE_HEIGHT;y++)
    {
        for(int x=0;x<IMAGE_WIDTH;x++)
        {
            output_image[y][x] = mask_image[y][x];
        }
    }

    /*
        Nur unbekannte Pixel ersetzen
    */
    for(int y=0;y<IMAGE_HEIGHT;y++)
    {
        for(int x=0;x<IMAGE_WIDTH;x++)
        {
            if(output_image[y][x] != UNKNOWN_PIXEL)
                continue;

            float value = 0.0f;

            for(int k=0;k<MAX_FREQS;k++)
            {
                float angle =
                    2.0f * PI *
                    (
                        ((float)components[k].u * x / IMAGE_WIDTH) +
                        ((float)components[k].v * y / IMAGE_HEIGHT)
                    )
                    + components[k].phase;

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

    fse_reconstruct_and_fill();

    /*
        Ausgabe kompletter Bildmatrix
    */

    for(int y=0;y<IMAGE_HEIGHT;y++)
    {
        for(int x=0;x<IMAGE_WIDTH;x++)
        {
            printf("%i\n", output_image[y][x]);
        }
        
    }

    return 0;
}