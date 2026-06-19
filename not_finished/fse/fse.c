#include <stdint.h>
#include <math.h>

#define WIDTH       64
#define HEIGHT      64

#define MAX_FREQS   32
#define LUT_SIZE    1024

#define PI 3.14159265358979323846f

typedef struct
{
    int16_t u;
    int16_t v;

    float amplitude;
    float phase;

} FSE_Component;

static float sin_lut[LUT_SIZE];
static float cos_lut[LUT_SIZE];

static FSE_Component components[MAX_FREQS];

static float residual[HEIGHT][WIDTH];
static float reconstruction[HEIGHT][WIDTH];

static inline int lut_index(float angle)
{
    while(angle < 0.0f)
        angle += 2.0f * PI;

    while(angle >= 2.0f * PI)
        angle -= 2.0f * PI;

    return (int)(angle * LUT_SIZE / (2.0f * PI));
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
        float a = 2.0f * PI * i / LUT_SIZE;

        sin_lut[i] = sinf(a);
        cos_lut[i] = cosf(a);
    }
}

static void copy_image( const float src[HEIGHT][WIDTH], float dst[HEIGHT][WIDTH]){
    for(int y=0;y<HEIGHT;y++)
    {
        for(int x=0;x<WIDTH;x++)
        {
            dst[y][x] = src[y][x];
        }
    }
}

static void clear_image(float img[HEIGHT][WIDTH]){
    for(int y=0;y<HEIGHT;y++)
    {
        for(int x=0;x<WIDTH;x++)
        {
            img[y][x] = 0.0f;
        }
    }
}

static void find_dominant_frequency( float img[HEIGHT][WIDTH], FSE_Component* comp){
    float best_energy = -1.0f;

    int best_u = 0;
    int best_v = 0;

    float best_re = 0.0f;
    float best_im = 0.0f;

    for(int v=0; v<HEIGHT/2; v++)
    {
        for(int u=0; u<WIDTH/2; u++)
        {
            float re = 0.0f;
            float im = 0.0f;

            for(int y=0;y<HEIGHT;y++)
            {
                for(int x=0;x<WIDTH;x++)
                {
                    float angle =
                        2.0f * PI *
                        (
                            (float)u*x/WIDTH +
                            (float)v*y/HEIGHT
                        );

                    re += img[y][x] * fast_cos(angle);
                    im -= img[y][x] * fast_sin(angle);
                }
            }

            float energy = re*re + im*im;

            if(energy > best_energy)
            {
                best_energy = energy;

                best_u = u;
                best_v = v;

                best_re = re;
                best_im = im;
            }
        }
    }

    comp->u = best_u;
    comp->v = best_v;

    comp->phase = atan2f(best_im, best_re);

    comp->amplitude =
        2.0f *
        sqrtf(best_energy) /
        (WIDTH * HEIGHT);
}

static void subtract_component(float residual_img[HEIGHT][WIDTH], const FSE_Component* c){
    for(int y=0;y<HEIGHT;y++)
    {
        for(int x=0;x<WIDTH;x++)
        {
            float angle =
                2.0f * PI *
                (
                    (float)c->u*x/WIDTH +
                    (float)c->v*y/HEIGHT
                ) +
                c->phase;

            residual_img[y][x] -=
                c->amplitude *
                fast_cos(angle);
        }
    }
}

static void add_component(float image[HEIGHT][WIDTH],const FSE_Component* c){
    for(int y=0;y<HEIGHT;y++)
    {
        for(int x=0;x<WIDTH;x++)
        {
            float angle =
                2.0f * PI *
                (
                    (float)c->u*x/WIDTH +
                    (float)c->v*y/HEIGHT
                ) +
                c->phase;

            image[y][x] +=
                c->amplitude *
                fast_cos(angle);
        }
    }
}

void fse_analyse(const float image[HEIGHT][WIDTH],int component_count){
    copy_image(image, residual);

    if(component_count > MAX_FREQS)
        component_count = MAX_FREQS;

    for(int i=0;i<component_count;i++)
    {
        find_dominant_frequency(
            residual,
            &components[i]);

        subtract_component(
            residual,
            &components[i]);
    }
}

void fse_reconstruct(float out[HEIGHT][WIDTH],int component_count){
    clear_image(out);

    if(component_count > MAX_FREQS)
        component_count = MAX_FREQS;

    for(int i=0;i<component_count;i++)
    {
        add_component(
            out,
            &components[i]);
    }
}

/*
Beispiel:

float image[HEIGHT][WIDTH];

fse_init();

fse_analyse(image, 32);

fse_reconstruct(
    reconstruction,
    32);

*/
