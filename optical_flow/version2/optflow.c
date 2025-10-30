// #include <stdio.h>
// #include <stdlib.h>
// #include <stdint.h>
// #include <math.h>
// #include <string.h>

// #define ALPHA 1.0f
// #define ITERATIONS 100

// // Horn–Schunck Optical Flow
// void horn_schunck(const float *I1, const float *I2, float *u, float *v, int w, int h, float alpha, int num_iter) {
//     float *Ix = malloc(sizeof(float)*w*h);
//     float *Iy = malloc(sizeof(float)*w*h);
//     float *It = malloc(sizeof(float)*w*h);
//     float *u_avg = malloc(sizeof(float)*w*h);
//     float *v_avg = malloc(sizeof(float)*w*h);

//     for (int y=1; y<h-1; y++) {
//         for (int x=1; x<w-1; x++) {
//             int i = y*w + x;
//             Ix[i] = 0.25f*((I1[i+1]-I1[i-1]) + (I2[i+1]-I2[i-1]));
//             Iy[i] = 0.25f*((I1[i+w]-I1[i-w]) + (I2[i+w]-I2[i-w]));
//             It[i] = I2[i] - I1[i];
//         }
//     }

//     for (int iter=0; iter<num_iter; iter++) {
//         for (int y=1; y<h-1; y++) {
//             for (int x=1; x<w-1; x++) {
//                 int i = y*w + x;
//                 u_avg[i] = 0.25f*(u[i-1]+u[i+1]+u[i-w]+u[i+w]);
//                 v_avg[i] = 0.25f*(v[i-1]+v[i+1]+v[i-w]+v[i+w]);
//                 float num = Ix[i]*u_avg[i] + Iy[i]*v_avg[i] + It[i];
//                 float den = alpha*alpha + Ix[i]*Ix[i] + Iy[i]*Iy[i];
//                 u[i] = u_avg[i] - Ix[i]*num/den;
//                 v[i] = v_avg[i] - Iy[i]*num/den;
//             }
//         }
//     }

//     free(Ix); free(Iy); free(It); free(u_avg); free(v_avg);
// }

// // Lade Video-Rohdaten
// uint8_t* load_raw_video(const char *filename, int *w, int *h, int *nframes) {
//     FILE *f = fopen(filename, "rb");
//     if (!f) { perror("fopen"); exit(1); }
//     fread(w, sizeof(uint32_t), 1, f);
//     fread(h, sizeof(uint32_t), 1, f);
//     fread(nframes, sizeof(uint32_t), 1, f);
//     size_t framesize = (*w)*(*h)*(*nframes);
//     uint8_t *frames = malloc(framesize);
//     fread(frames, 1, framesize, f);
//     fclose(f);
//     return frames;
// }

// // Speichere Flow in .raw-Datei
// void save_flow_raw(const char *filename, float *u, float *v, int w, int h, int nframes) {
//     FILE *f = fopen(filename, "wb");
//     if (!f) { perror("fopen"); exit(1); }
//     fwrite(&w, sizeof(uint32_t), 1, f);
//     fwrite(&h, sizeof(uint32_t), 1, f);
//     fwrite(&nframes, sizeof(uint32_t), 1, f);
//     fwrite(u, sizeof(float), w*h*nframes, f);
//     fwrite(v, sizeof(float), w*h*nframes, f);
//     fclose(f);
// }

// int main(int argc, char **argv) {
//     if (argc < 3) {
//         printf("Verwendung: %s <input.raw> <output.raw>\n", argv[0]);
//         return 1;
//     }

//     int w,h,nframes;
//     uint8_t *frames = load_raw_video(argv[1], &w, &h, &nframes);
//     printf("Video geladen: %dx%d, %d Frames\n", w, h, nframes);

//     float *I1 = malloc(sizeof(float)*w*h);
//     float *I2 = malloc(sizeof(float)*w*h);
//     float *u = malloc(sizeof(float)*w*h*(nframes-1));
//     float *v = malloc(sizeof(float)*w*h*(nframes-1));

//     for (int f=0; f<nframes-1; f++) {
//         for (int i=0; i<w*h; i++) {
//             I1[i] = frames[f*w*h + i] / 255.0f;
//             I2[i] = frames[(f+1)*w*h + i] / 255.0f;
//         }

//         float *uf = u + f*w*h;
//         float *vf = v + f*w*h;
//         memset(uf, 0, sizeof(float)*w*h);
//         memset(vf, 0, sizeof(float)*w*h);

//         horn_schunck(I1, I2, uf, vf, w, h, ALPHA, ITERATIONS);
//         printf("Frame %d / %d berechnet\n", f+1, nframes-1);
//     }

//     save_flow_raw(argv[2], u, v, w, h, nframes-1);
//     printf("Flow gespeichert: %s\n", argv[2]);

//     free(frames); free(I1); free(I2); free(u); free(v);
//     return 0;
// }


/*
 lucas_kanade_from_raw.c

 Dense Lucas-Kanade optical flow (windowed) on RAW video frames.

 Usage:
   ./lucas_kanade_from_raw input_frames.raw flow_output.raw [window_radius]

 Input .raw format:
   uint32_t width
   uint32_t height
   uint32_t frame_count
   uint8_t frame1[w*h]
   uint8_t frame2[w*h]
   ...

 Output .raw format:
   uint32_t width
   uint32_t height
   uint32_t frame_count  (this equals input_frame_count - 1)
   float  u_frame1[w*h]
   float  u_frame2[w*h]
   ...
   float  v_frame1[w*h]
   float  v_frame2[w*h]
   ...

 Notes:
 - Uses integral images for fast box-sums over window (O(1) per pixel).
 - Default window_radius = 2 (i.e. 5x5 window). Increase for more robustness.
 - Uses floats. For RISC-V without FPU, compile with -msoft-float or consider fixed-point.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#ifndef CLAMP
#define CLAMP(x,a,b) ((x)<(a)?(a):((x)>(b)?(b):(x)))
#endif

// read raw frames: returns pointer to allocated frames (uint8_t array of size w*h*nframes)
uint8_t *load_raw_video(const char *fname, int *w, int *h, int *nframes) {
    FILE *f = fopen(fname, "rb");
    if(!f){ perror("fopen input"); exit(1); }
    uint32_t wu, hu, nu;
    if(fread(&wu, sizeof(uint32_t), 1, f) != 1) { perror("read w"); exit(1); }
    if(fread(&hu, sizeof(uint32_t), 1, f) != 1) { perror("read h"); exit(1); }
    if(fread(&nu, sizeof(uint32_t), 1, f) != 1) { perror("read n"); exit(1); }
    *w = (int)wu; *h = (int)hu; *nframes = (int)nu;
    size_t framesize = (size_t)(*w)*(*h)*(*nframes);
    uint8_t *buf = (uint8_t*)malloc(framesize);
    if(!buf){ fprintf(stderr,"Out of memory\n"); exit(1); }
    size_t got = fread(buf, 1, framesize, f);
    if(got != framesize){ fprintf(stderr, "Unexpected input size: got %zu expected %zu\n", got, framesize); exit(1); }
    fclose(f);
    return buf;
}

void save_flow_raw(const char *fname, float *u_all, float *v_all, int w, int h, int nframes) {
    FILE *f = fopen(fname, "wb");
    if(!f){ perror("fopen output"); exit(1); }
    uint32_t wu = (uint32_t)w, hu = (uint32_t)h, nu = (uint32_t)nframes;
    fwrite(&wu, sizeof(uint32_t), 1, f);
    fwrite(&hu, sizeof(uint32_t), 1, f);
    fwrite(&nu, sizeof(uint32_t), 1, f);
    // write all u then all v (as in design)
    size_t total = (size_t)w*h*nframes;
    fwrite(u_all, sizeof(float), total, f);
    fwrite(v_all, sizeof(float), total, f);
    fclose(f);
}

// compute gradients Ix, Iy, It (float arrays) for image pair I1,I2 (both float, normalized [0..1])
void compute_gradients(const float *I1, const float *I2, float *Ix, float *Iy, float *It, int w, int h) {
    // central differences for interior, forward/backward for borders
    for(int y=0;y<h;y++){
        for(int x=0;x<w;x++){
            int i = y*w + x;
            int xl = (x==0) ? x : x-1;
            int xr = (x==w-1) ? x : x+1;
            int yu = (y==0) ? y : y-1;
            int yd = (y==h-1) ? y : y+1;
            Ix[i] = 0.5f * (I1[y*w + xr] - I1[y*w + xl] + I2[y*w + xr] - I2[y*w + xl]);
            Iy[i] = 0.5f * (I1[yd*w + x] - I1[yu*w + x] + I2[yd*w + x] - I2[yu*w + x]);
            It[i] = I2[i] - I1[i];
        }
    }
}

// build integral image for array A (size w*h)
// integral has size (h+1)*(w+1), row-major, with first row/col zero
void build_integral(const float *A, float *intg, int w, int h) {
    int gw = w+1;
    // zero first row
    for(int x=0;x<gw;x++) intg[x] = 0.0f;
    for(int y=1;y<=h;y++){
        float rowsum = 0.0f;
        int row_idx = y*gw;
        int src_row = (y-1)*w;
        intg[row_idx + 0] = 0.0f;
        for(int x=1;x<=w;x++){
            rowsum += A[src_row + (x-1)];
            intg[row_idx + x] = intg[row_idx - gw + x] + rowsum;
        }
    }
}

// sum over box [x0..x1] [y0..y1] inclusive using integral image
static inline float integral_sum(const float *intg, int gw, int x0, int y0, int x1, int y1) {
    // coordinates in integral: +1 offset
    int Ax = x0, Ay = y0, Bx = x1+1, By = y1+1;
    // intg[row*gw + col]
    float A = intg[Ay*gw + Ax];
    float B = intg[Ay*gw + Bx];
    float C = intg[By*gw + Ax];
    float D = intg[By*gw + Bx];
    return D - B - C + A;
}

int main(int argc, char **argv) {
    if(argc < 3) {
        fprintf(stderr, "Usage: %s input_frames.raw flow_output.raw [window_radius]\n", argv[0]);
        return 1;
    }
    const char *infile = argv[1];
    const char *outfile = argv[2];
    int radius = 2; // default window radius -> window size = 5
    if(argc >= 4) radius = atoi(argv[3]);
    if(radius < 1) radius = 1;

    int w,h,nframes;
    uint8_t *raw = load_raw_video(infile, &w, &h, &nframes);
    if(nframes < 2) { fprintf(stderr,"Need >=2 frames\n"); return 1; }
    printf("Loaded: %d frames, size %dx%d, window radius %d\n", nframes, w, h, radius);

    const int N = w*h;
    // buffers per frame-pair
    float *I1 = (float*)malloc(sizeof(float)*N);
    float *I2 = (float*)malloc(sizeof(float)*N);
    float *Ix = (float*)malloc(sizeof(float)*N);
    float *Iy = (float*)malloc(sizeof(float)*N);
    float *It = (float*)malloc(sizeof(float)*N);

    // product arrays for integral images
    float *Ix2 = (float*)malloc(sizeof(float)*N);
    float *Iy2 = (float*)malloc(sizeof(float)*N);
    float *Ixy = (float*)malloc(sizeof(float)*N);
    float *IxIt = (float*)malloc(sizeof(float)*N);
    float *IyIt = (float*)malloc(sizeof(float)*N);

    // integral images (dimensions (h+1)*(w+1))
    int gw = w + 1;
    int gh = h + 1;
    size_t gint_size = (size_t)gw * gh;
    float *int_Ix2 = (float*)malloc(sizeof(float)*gint_size);
    float *int_Iy2 = (float*)malloc(sizeof(float)*gint_size);
    float *int_Ixy = (float*)malloc(sizeof(float)*gint_size);
    float *int_IxIt = (float*)malloc(sizeof(float)*gint_size);
    float *int_IyIt = (float*)malloc(sizeof(float)*gint_size);

    if(!I1||!I2||!Ix||!Iy||!It||!Ix2||!Iy2||!Ixy||!IxIt||!IyIt||!int_Ix2) {
        fprintf(stderr,"Malloc failed\n"); return 1;
    }

    int out_frames = nframes - 1;
    float *u_all = (float*)malloc(sizeof(float)*N*out_frames);
    float *v_all = (float*)malloc(sizeof(float)*N*out_frames);
    if(!u_all || !v_all){ fprintf(stderr,"Out of memory alloc out\n"); return 1; }

    // process each pair
    for(int f=0; f<out_frames; ++f){
        // convert to float normalized
        uint8_t *f1 = raw + (size_t)f * N;
        uint8_t *f2 = raw + (size_t)(f+1) * N;
        for(int i=0;i<N;i++){
            I1[i] = (float)f1[i] / 255.0f;
            I2[i] = (float)f2[i] / 255.0f;
        }

        // gradients
        compute_gradients(I1, I2, Ix, Iy, It, w, h);

        // products
        for(int i=0;i<N;i++){
            float ix = Ix[i], iy = Iy[i], it = It[i];
            Ix2[i] = ix * ix;
            Iy2[i] = iy * iy;
            Ixy[i] = ix * iy;
            IxIt[i] = ix * it;
            IyIt[i] = iy * it;
        }

        // build integrals
        build_integral(Ix2, int_Ix2, w, h);
        build_integral(Iy2, int_Iy2, w, h);
        build_integral(Ixy, int_Ixy, w, h);
        build_integral(IxIt, int_IxIt, w, h);
        build_integral(IyIt, int_IyIt, w, h);

        // for each pixel compute sums over window and solve
        int winsize_x0 = radius, winsize_y0 = radius;
        for(int y=0;y<h;y++){
            int y0 = CLAMP(y - radius, 0, h-1);
            int y1 = CLAMP(y + radius, 0, h-1);
            for(int x=0;x<w;x++){
                int x0 = CLAMP(x - radius, 0, w-1);
                int x1 = CLAMP(x + radius, 0, w-1);
                // sums using integral images; integral gw = w+1
                float S_Ix2 = integral_sum(int_Ix2, gw, x0, y0, x1, y1);
                float S_Iy2 = integral_sum(int_Iy2, gw, x0, y0, x1, y1);
                float S_Ixy = integral_sum(int_Ixy, gw, x0, y0, x1, y1);
                float S_IxIt = integral_sum(int_IxIt, gw, x0, y0, x1, y1);
                float S_IyIt = integral_sum(int_IyIt, gw, x0, y0, x1, y1);

                // Solve: [Sxx Sxy; Sxy Syy] * [u; v] = -[SxIt; SyIt]
                float det = S_Ix2 * S_Iy2 - S_Ixy * S_Ixy;
                float uval = 0.0f, vval = 0.0f;
                const float eps = 1e-6f;
                if(fabsf(det) > eps) {
                    float inv_det = 1.0f / det;
                    uval = (-S_Iy2 * S_IxIt + S_Ixy * S_IyIt) * inv_det;
                    vval = ( S_Ixy * S_IxIt - S_Ix2 * S_IyIt) * inv_det;
                } else {
                    uval = 0.0f; vval = 0.0f;
                }
                int idx = y*w + x;
                u_all[(size_t)f * N + idx] = uval;
                v_all[(size_t)f * N + idx] = vval;
            }
        }
        printf("Processed frame %d/%d\n", f+1, out_frames);
    }

    // save
    save_flow_raw(outfile, u_all, v_all, w, h, out_frames);
    printf("Saved flow to %s\n", outfile);

    // free
    free(raw);
    free(I1); free(I2); free(Ix); free(Iy); free(It);
    free(Ix2); free(Iy2); free(Ixy); free(IxIt); free(IyIt);
    free(int_Ix2); free(int_Iy2); free(int_Ixy); free(int_IxIt); free(int_IyIt);
    free(u_all); free(v_all);
    return 0;
}
