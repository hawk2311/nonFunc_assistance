// // #include <stdio.h>
// // #include <stdlib.h>
// // #include <stdint.h>
// // #include <math.h>
// // #include <string.h>

// // #define ALPHA 1.0f
// // #define ITERATIONS 100

// // // Horn–Schunck Optical Flow
// // void horn_schunck(const float *I1, const float *I2, float *u, float *v, int w, int h, float alpha, int num_iter) {
// //     float *Ix = malloc(sizeof(float)*w*h);
// //     float *Iy = malloc(sizeof(float)*w*h);
// //     float *It = malloc(sizeof(float)*w*h);
// //     float *u_avg = malloc(sizeof(float)*w*h);
// //     float *v_avg = malloc(sizeof(float)*w*h);

// //     for (int y=1; y<h-1; y++) {
// //         for (int x=1; x<w-1; x++) {
// //             int i = y*w + x;
// //             Ix[i] = 0.25f*((I1[i+1]-I1[i-1]) + (I2[i+1]-I2[i-1]));
// //             Iy[i] = 0.25f*((I1[i+w]-I1[i-w]) + (I2[i+w]-I2[i-w]));
// //             It[i] = I2[i] - I1[i];
// //         }
// //     }

// //     for (int iter=0; iter<num_iter; iter++) {
// //         for (int y=1; y<h-1; y++) {
// //             for (int x=1; x<w-1; x++) {
// //                 int i = y*w + x;
// //                 u_avg[i] = 0.25f*(u[i-1]+u[i+1]+u[i-w]+u[i+w]);
// //                 v_avg[i] = 0.25f*(v[i-1]+v[i+1]+v[i-w]+v[i+w]);
// //                 float num = Ix[i]*u_avg[i] + Iy[i]*v_avg[i] + It[i];
// //                 float den = alpha*alpha + Ix[i]*Ix[i] + Iy[i]*Iy[i];
// //                 u[i] = u_avg[i] - Ix[i]*num/den;
// //                 v[i] = v_avg[i] - Iy[i]*num/den;
// //             }
// //         }
// //     }

// //     free(Ix); free(Iy); free(It); free(u_avg); free(v_avg);
// // }

// // // Lade Video-Rohdaten
// // uint8_t* load_raw_video(const char *filename, int *w, int *h, int *nframes) {
// //     FILE *f = fopen(filename, "rb");
// //     if (!f) { perror("fopen"); exit(1); }
// //     fread(w, sizeof(uint32_t), 1, f);
// //     fread(h, sizeof(uint32_t), 1, f);
// //     fread(nframes, sizeof(uint32_t), 1, f);
// //     size_t framesize = (*w)*(*h)*(*nframes);
// //     uint8_t *frames = malloc(framesize);
// //     fread(frames, 1, framesize, f);
// //     fclose(f);
// //     return frames;
// // }

// // // Speichere Flow in .raw-Datei
// // void save_flow_raw(const char *filename, float *u, float *v, int w, int h, int nframes) {
// //     FILE *f = fopen(filename, "wb");
// //     if (!f) { perror("fopen"); exit(1); }
// //     fwrite(&w, sizeof(uint32_t), 1, f);
// //     fwrite(&h, sizeof(uint32_t), 1, f);
// //     fwrite(&nframes, sizeof(uint32_t), 1, f);
// //     fwrite(u, sizeof(float), w*h*nframes, f);
// //     fwrite(v, sizeof(float), w*h*nframes, f);
// //     fclose(f);
// // }

// // int main(int argc, char **argv) {
// //     if (argc < 3) {
// //         printf("Verwendung: %s <input.raw> <output.raw>\n", argv[0]);
// //         return 1;
// //     }

// //     int w,h,nframes;
// //     uint8_t *frames = load_raw_video(argv[1], &w, &h, &nframes);
// //     printf("Video geladen: %dx%d, %d Frames\n", w, h, nframes);

// //     float *I1 = malloc(sizeof(float)*w*h);
// //     float *I2 = malloc(sizeof(float)*w*h);
// //     float *u = malloc(sizeof(float)*w*h*(nframes-1));
// //     float *v = malloc(sizeof(float)*w*h*(nframes-1));

// //     for (int f=0; f<nframes-1; f++) {
// //         for (int i=0; i<w*h; i++) {
// //             I1[i] = frames[f*w*h + i] / 255.0f;
// //             I2[i] = frames[(f+1)*w*h + i] / 255.0f;
// //         }

// //         float *uf = u + f*w*h;
// //         float *vf = v + f*w*h;
// //         memset(uf, 0, sizeof(float)*w*h);
// //         memset(vf, 0, sizeof(float)*w*h);

// //         horn_schunck(I1, I2, uf, vf, w, h, ALPHA, ITERATIONS);
// //         printf("Frame %d / %d berechnet\n", f+1, nframes-1);
// //     }

// //     save_flow_raw(argv[2], u, v, w, h, nframes-1);
// //     printf("Flow gespeichert: %s\n", argv[2]);

// //     free(frames); free(I1); free(I2); free(u); free(v);
// //     return 0;
// // }


// /*
//  lucas_kanade_from_raw.c

//  Dense Lucas-Kanade optical flow (windowed) on RAW video frames.

//  Usage:
//    ./lucas_kanade_from_raw input_frames.raw flow_output.raw [window_radius]

//  Input .raw format:
//    uint32_t width
//    uint32_t height
//    uint32_t frame_count
//    uint8_t frame1[w*h]
//    uint8_t frame2[w*h]
//    ...

//  Output .raw format:
//    uint32_t width
//    uint32_t height
//    uint32_t frame_count  (this equals input_frame_count - 1)
//    float  u_frame1[w*h]
//    float  u_frame2[w*h]
//    ...
//    float  v_frame1[w*h]
//    float  v_frame2[w*h]
//    ...

//  Notes:
//  - Uses integral images for fast box-sums over window (O(1) per pixel).
//  - Default window_radius = 2 (i.e. 5x5 window). Increase for more robustness.
//  - Uses floats. For RISC-V without FPU, compile with -msoft-float or consider fixed-point.
// */

// /*
//  lucas_kanade_baremetal.c
//  Dense Lucas-Kanade for RAW frames.
//  Two modes:
//    - HOST mode (compile with -DHOST): reads/writes .raw via fopen (for test on host or RISC-V Linux)
//    - Bare-metal mode (default): include "video_data.h" (xxd -i ...) which defines `unsigned char video_data[]` and `unsigned int video_data_len`;
//        program reads header+frames from that array and emits the flow_output.raw as binary bytes via SBI console_putchar (bytewise).
//  Usage (HOST):
//    gcc -O3 -DHOST -o lk_host lucas_kanade_baremetal.c -lm
//    ./lk_host video_frames.raw flow_output.raw [window_radius]
//  Usage (Bare-metal):
//    # create video_data.h from video_frames.raw:
//    xxd -i video_frames.raw > video_data.h
//    # compile:
//    riscv64-unknown-elf-gcc -O3 -march=rv64imac -mabi=lp64 -o lk_bare.elf lucas_kanade_baremetal.c
//    # run under spike and capture output:
//    spike lk_bare.elf > flow_output.raw
//  Note: If using pk (proxy kernel) on spike, you can also use `spike pk lk_bare.elf > flow_output.raw`
 
//  Output .raw format:
//    uint32_t width, uint32_t height, uint32_t frame_count (nframes_out)
//    float u_frame0[w*h] ... float u_frameN[w*h]
//    float v_frame0[w*h] ... float v_frameN[w*h]
// */

// /*
//  lucas_kanade_baremetal.c
//  Dense Lucas-Kanade for RAW frames.
//  Two modes:
//    - HOST mode (compile with -DHOST): reads/writes .raw via fopen (for test on host or RISC-V Linux)
//    - Bare-metal mode (default): include "video_data.h" (xxd -i ...) which defines `unsigned char video_data[]` and `unsigned int video_data_len`;
//        program reads header+frames from that array and emits the flow_output.raw as binary bytes via SBI console_putchar (bytewise).
//  Usage (HOST):
//    gcc -O3 -DHOST -o lk_host lucas_kanade_baremetal.c -lm
//    ./lk_host video_frames.raw flow_output.raw [window_radius]
//  Usage (Bare-metal):
//    # create video_data.h from video_frames.raw:
//    xxd -i video_frames.raw > video_data.h
//    # compile:
//    riscv64-unknown-elf-gcc -O3 -march=rv64imac -mabi=lp64 -o lk_bare.elf lucas_kanade_baremetal.c
//    # run under spike and capture output:
//    spike lk_bare.elf > flow_output.raw
//  Note: If using pk (proxy kernel) on spike, you can also use `spike pk lk_bare.elf > flow_output.raw`
 
//  Output .raw format:
//    uint32_t width, uint32_t height, uint32_t frame_count (nframes_out)
//    float u_frame0[w*h] ... float u_frameN[w*h]
//    float v_frame0[w*h] ... float v_frameN[w*h]
// */

// #include <stdint.h>
// #include <stddef.h>

// #ifdef HOST
//   #include <stdio.h>
//   #include <stdlib.h>
//   #include <string.h>
//   #include <math.h>
// #else
//   /* Bare-metal: provide ecall console putchar (SBI) */
//   static inline void sbicon_putc(unsigned char c){
//       register long a0 asm("a0") = (long)c;
//       register long a7 asm("a7") = 1; /* sbi_console_putchar legacy */
//       asm volatile("ecall" : "+r"(a0) : "r"(a7) : "memory");
//   }
// #endif

// #ifndef HOST
//   /* user must generate video_data.h with: xxd -i video_frames.raw > video_data.h
//      it must define: unsigned char video_data[]; unsigned int video_data_len; */
//   #include "video_data.h"
// #endif

// /* configurable maximums for bare-metal static allocation */
// #ifndef MAX_W
//   #define MAX_W 320
// #endif
// #ifndef MAX_H
//   #define MAX_H 240
// #endif
// #ifndef MAX_FRAMES
//   #define MAX_FRAMES 200
// #endif

// /* helper for little-endian read/write */
// static inline uint32_t read_u32_le_from_buf(const unsigned char *b){
//     return (uint32_t)b[0] | ((uint32_t)b[1]<<8) | ((uint32_t)b[2]<<16) | ((uint32_t)b[3]<<24);
// }
// static inline void write_u32_le_to_buf(unsigned char *b, uint32_t v){
//     b[0] = v & 0xFF; b[1] = (v>>8)&0xFF; b[2] = (v>>16)&0xFF; b[3] = (v>>24)&0xFF;
// }
// static inline void write_float_le_bytes(unsigned char *b, float f){
//     union { float f; uint32_t u; } x; x.f = f;
//     write_u32_le_to_buf(b, x.u);
// }

// /* Integral image helpers (float) */
// /* integral image dims = (h+1) x (w+1) */
// static inline void build_integral_f(const float *src, float *intg, int w, int h){
//     int gw = w+1;
//     // first row zero
//     for(int x=0;x<gw;x++) intg[x] = 0.0f;
//     for(int y=1;y<=h;y++){
//         float rsum = 0.0f;
//         int dst_row = y*gw;
//         int src_row = (y-1)*w;
//         intg[dst_row + 0] = 0.0f;
//         for(int x=1;x<=w;x++){
//             rsum += src[src_row + (x-1)];
//             intg[dst_row + x] = intg[dst_row - gw + x] + rsum;
//         }
//     }
// }
// static inline float integral_sum_f(const float *intg, int gw, int x0, int y0, int x1, int y1){
//     // coords in intg: +1 shift
//     int Ax = x0;
//     int Ay = y0;
//     int Bx = x1 + 1;
//     int By = y1 + 1;
//     float A = intg[Ay*gw + Ax];
//     float B = intg[Ay*gw + Bx];
//     float C = intg[By*gw + Ax];
//     float D = intg[By*gw + Bx];
//     return D - B - C + A;
// }

// /* ---------- MAIN PROCESSING LOGIC ---------- */

// /* We'll allocate static buffers (bare-metal friendly). For HOST mode we may malloc. */
// static float I1_buf[MAX_W * MAX_H];
// static float I2_buf[MAX_W * MAX_H];
// static float Ix_buf[MAX_W * MAX_H];
// static float Iy_buf[MAX_W * MAX_H];
// static float It_buf[MAX_W * MAX_H];
// static float Ix2_buf[MAX_W * MAX_H];
// static float Iy2_buf[MAX_W * MAX_H];
// static float Ixy_buf[MAX_W * MAX_H];
// static float IxIt_buf[MAX_W * MAX_H];
// static float IyIt_buf[MAX_W * MAX_H];

// /* integral images size (w+1)*(h+1) - we allocate max */
// static float int_Ix2[(MAX_W+1)*(MAX_H+1)];
// static float int_Iy2[(MAX_W+1)*(MAX_H+1)];
// static float int_Ixy[(MAX_W+1)*(MAX_H+1)];
// static float int_IxIt[(MAX_W+1)*(MAX_H+1)];
// static float int_IyIt[(MAX_W+1)*(MAX_H+1)];

// /* small utility to send bytes in bare-metal (SBI) */
// #ifndef HOST
// static void send_bytes_bare(const unsigned char *buf, size_t len){
//     for(size_t i=0;i<len;i++) sbicon_putc(buf[i]);
// }
// #endif

// #ifdef HOST
// /* Host helpers (file IO) */
// #include <stdlib.h>
// #include <stdio.h>

// static uint8_t *load_raw_host(const char *fname, int *w, int *h, int *nframes){
//     FILE *f = fopen(fname,"rb");
//     if(!f){ perror("open"); return NULL; }
//     uint32_t wu, hu, nu;
//     if(fread(&wu,4,1,f)!=1){ fclose(f); return NULL; }
//     if(fread(&hu,4,1,f)!=1){ fclose(f); return NULL; }
//     if(fread(&nu,4,1,f)!=1){ fclose(f); return NULL; }
//     *w = (int)wu; *h = (int)hu; *nframes = (int)nu;
//     size_t framesz = (size_t)(*w)*(*h)*(*nframes);
//     uint8_t *buf = malloc(framesz);
//     if(!buf){ fclose(f); return NULL; }
//     if(fread(buf,1,framesz,f)!=framesz){ free(buf); fclose(f); return NULL; }
//     fclose(f);
//     return buf;
// }

// static int save_flow_host(const char *fname, float *u_all, float *v_all, int w, int h, int nframes){
//     FILE *f = fopen(fname,"wb");
//     if(!f){ perror("open out"); return -1; }
//     uint32_t wu = (uint32_t)w, hu = (uint32_t)h, nu = (uint32_t)nframes;
//     fwrite(&wu,4,1,f); fwrite(&hu,4,1,f); fwrite(&nu,4,1,f);
//     size_t Nall = (size_t)w*h*nframes;
//     fwrite(u_all, sizeof(float), Nall, f);
//     fwrite(v_all, sizeof(float), Nall, f);
//     fclose(f);
//     return 0;
// }
// #endif

// /* compute gradients (central differences) */
// static void compute_gradients(const float *I1, const float *I2, float *Ix, float *Iy, float *It, int w, int h){
//     for(int y=0;y<h;y++){
//         for(int x=0;x<w;x++){
//             int xl = (x==0)?x:x-1;
//             int xr = (x==w-1)?x:x+1;
//             int yu = (y==0)?y:y-1;
//             int yd = (y==h-1)?y:y+1;
//             int idx = y*w + x;
//             Ix[idx] = 0.5f * (I1[y*w + xr] - I1[y*w + xl] + I2[y*w + xr] - I2[y*w + xl]);
//             Iy[idx] = 0.5f * (I1[yd*w + x] - I1[yu*w + x] + I2[yd*w + x] - I2[yu*w + x]);
//             It[idx] = I2[idx] - I1[idx];
//         }
//     }
// }

// #ifdef HOST
// /* HOST main path */
// int main(int argc, char **argv){
//     if(argc < 3){
//         fprintf(stderr,"Usage: %s input_frames.raw flow_output.raw [window_radius]\n", argv[0]);
//         return 1;
//     }
//     const char *infile = argv[1];
//     const char *outfile = argv[2];
//     int radius = 2;
//     if(argc>=4) radius = atoi(argv[3]);
//     if(radius < 1) radius = 1;

//     int w,h,nframes;
//     uint8_t *raw = load_raw_host(infile,&w,&h,&nframes);
//     if(!raw){ fprintf(stderr,"Failed to load input\n"); return 1; }
//     int N = w*h;
//     int out_frames = nframes - 1;
//     if(out_frames <= 0){ fprintf(stderr,"Need >=2 frames\n"); return 1; }

//     float *u_all = (float*)malloc(sizeof(float)*N*out_frames);
//     float *v_all = (float*)malloc(sizeof(float)*N*out_frames);
//     if(!u_all || !v_all){ fprintf(stderr,"Out of memory\n"); return 1; }

//     for(int f=0; f<out_frames; ++f){
//         /* fill I1,I2 */
//         uint8_t *p1 = raw + (size_t)f * N;
//         uint8_t *p2 = raw + (size_t)(f+1) * N;
//         for(int i=0;i<N;i++){ I1_buf[i] = (float)p1[i] / 255.0f; I2_buf[i] = (float)p2[i] / 255.0f; }

//         compute_gradients(I1_buf, I2_buf, Ix_buf, Iy_buf, It_buf, w, h);

//         /* compute products */
//         for(int i=0;i<N;i++){
//             float ix = Ix_buf[i], iy = Iy_buf[i], it = It_buf[i];
//             Ix2_buf[i] = ix*ix;
//             Iy2_buf[i] = iy*iy;
//             Ixy_buf[i] = ix*iy;
//             IxIt_buf[i] = ix*it;
//             IyIt_buf[i] = iy*it;
//         }
//         /* build integrals */
//         int gw = w+1;
//         build_integral_f(Ix2_buf, int_Ix2, w, h);
//         build_integral_f(Iy2_buf, int_Iy2, w, h);
//         build_integral_f(Ixy_buf, int_Ixy, w, h);
//         build_integral_f(IxIt_buf, int_IxIt, w, h);
//         build_integral_f(IyIt_buf, int_IyIt, w, h);

//         /* compute per-pixel */
//         for(int y=0;y<h;y++){
//             int y0 = (y - radius) < 0 ? 0 : (y - radius);
//             int y1 = (y + radius) >= h ? (h-1) : (y + radius);
//             for(int x=0;x<w;x++){
//                 int x0 = (x - radius) < 0 ? 0 : (x - radius);
//                 int x1 = (x + radius) >= w ? (w-1) : (x + radius);

//                 float S_Ix2 = integral_sum_f(int_Ix2, gw, x0, y0, x1, y1);
//                 float S_Iy2 = integral_sum_f(int_Iy2, gw, x0, y0, x1, y1);
//                 float S_Ixy = integral_sum_f(int_Ixy, gw, x0, y0, x1, y1);
//                 float S_IxIt = integral_sum_f(int_IxIt, gw, x0, y0, x1, y1);
//                 float S_IyIt = integral_sum_f(int_IyIt, gw, x0, y0, x1, y1);

//                 float det = S_Ix2 * S_Iy2 - S_Ixy * S_Ixy;
//                 float uval = 0.0f, vval = 0.0f;
//                 const float eps = 1e-6f;
//                 if(det > eps || det < -eps){
//                     float inv = 1.0f / det;
//                     uval = (-S_Iy2 * S_IxIt + S_Ixy * S_IyIt) * inv;
//                     vval = ( S_Ixy * S_IxIt - S_Ix2 * S_IyIt) * inv;
//                 }
//                 int idx = y*w + x;
//                 u_all[(size_t)f * N + idx] = uval;
//                 v_all[(size_t)f * N + idx] = vval;
//             }
//         }
//         fprintf(stderr,"Processed %d/%d\n", f+1, out_frames);
//     }

//     /* save */
//     if(save_flow_host(outfile, u_all, v_all, w, h, out_frames) != 0){
//         fprintf(stderr,"Failed save\n");
//     } else {
//         fprintf(stderr,"Saved %s\n", outfile);
//     }

//     free(raw);
//     free(u_all); free(v_all);
//     return 0;
// }
// #else
// /* ---------- BARE-METAL entry ---------- */
// /* read header and frames from embedded video_data[] and stream binary output via sbicon_putc */
// int main(void){
//     const unsigned char *p = video_data;
//     unsigned int len = video_data_len;
//     if(len < 12){
//         const char *err = "video_data too small\n"; for(const char *s=err; *s; ++s) sbicon_putc(*s);
//         while(1) asm volatile("wfi");
//     }
//     uint32_t w = read_u32_le_from_buf(p); p+=4;
//     uint32_t h = read_u32_le_from_buf(p); p+=4;
//     uint32_t nframes = read_u32_le_from_buf(p); p+=4;
//     const unsigned char *frames_ptr = p;
//     size_t expected = (size_t)w * h * nframes;
//     if(expected + 12 > len){
//         const char *err = "video_data length mismatch\n"; for(const char *s=err; *s; ++s) sbicon_putc(*s);
//         while(1) asm volatile("wfi");
//     }
//     int radius = 2; /* default; could be changed via compile-time macro */
//     int N = (int)w*(int)h;
//     int out_frames = (int)nframes - 1;
//     if(out_frames <= 0){
//         const char *err = "Need >=2 frames\n"; for(const char *s=err; *s; ++s) sbicon_putc(*s);
//         while(1) asm volatile("wfi");
//     }

//     /* static buffers used (I1/I2 etc above) */
//     for(int f=0; f<out_frames; ++f){
//         const unsigned char *f1 = frames_ptr + (size_t)f * N;
//         const unsigned char *f2 = frames_ptr + (size_t)(f+1) * N;
//         for(int i=0;i<N;i++){ I1_buf[i] = (float)f1[i] / 255.0f; I2_buf[i] = (float)f2[i] / 255.0f; }
//         compute_gradients(I1_buf, I2_buf, Ix_buf, Iy_buf, It_buf, (int)w, (int)h);
//         for(int i=0;i<N;i++){
//             float ix = Ix_buf[i], iy = Iy_buf[i], it = It_buf[i];
//             Ix2_buf[i] = ix*ix;
//             Iy2_buf[i] = iy*iy;
//             Ixy_buf[i] = ix*iy;
//             IxIt_buf[i] = ix*it;
//             IyIt_buf[i] = iy*it;
//         }
//         int gw = (int)w + 1;
//         /* build integrals into statically allocated arrays */
//         build_integral_f(Ix2_buf, int_Ix2, (int)w, (int)h);
//         build_integral_f(Iy2_buf, int_Iy2, (int)w, (int)h);
//         build_integral_f(Ixy_buf, int_Ixy, (int)w, (int)h);
//         build_integral_f(IxIt_buf, int_IxIt, (int)w, (int)h);
//         build_integral_f(IyIt_buf, int_IyIt, (int)w, (int)h);

//         /* For bare-metal we stream the results as we compute them to save RAM:
//            We'll first write header after processing all frames (but easier: collect in RAM is heavy).
//            To keep it simple, we'll stream header then (u,v) per frame.
//            We'll first send header: w,h,out_frames as 3 little-endian uint32.
//         */
//         /* We'll stream for each frame: first u-frame (w*h floats), then v-frame. */
//         /* To make the receiver compatible with save_flow_host format, we must send header FIRST.
//            So for bare-metal, send header before computation, then send per-frame data as floats.
//         */
//         if(f==0){
//             unsigned char hdr[12];
//             write_u32_le_to_buf(&hdr[0], (uint32_t)w);
//             write_u32_le_to_buf(&hdr[4], (uint32_t)h);
//             write_u32_le_to_buf(&hdr[8], (uint32_t)out_frames);
//             send_bytes_bare(hdr,12);
//         }

//         for(int y=0;y<(int)h;y++){
//             int y0 = (y - radius) < 0 ? 0 : (y - radius);
//             int y1 = (y + radius) >= (int)h ? (int)h-1 : (y + radius);
//             for(int x=0;x<(int)w;x++){
//                 int x0 = (x - radius) < 0 ? 0 : (x - radius);
//                 int x1 = (x + radius) >= (int)w ? (int)w-1 : (x + radius);
//                 float S_Ix2 = integral_sum_f(int_Ix2, gw, x0, y0, x1, y1);
//                 float S_Iy2 = integral_sum_f(int_Iy2, gw, x0, y0, x1, y1);
//                 float S_Ixy = integral_sum_f(int_Ixy, gw, x0, y0, x1, y1);
//                 float S_IxIt = integral_sum_f(int_IxIt, gw, x0, y0, x1, y1);
//                 float S_IyIt = integral_sum_f(int_IyIt, gw, x0, y0, x1, y1);
//                 float det = S_Ix2 * S_Iy2 - S_Ixy * S_Ixy;
//                 float uval=0.0f, vval=0.0f;
//                 const float eps = 1e-6f;
//                 if(det > eps || det < -eps){
//                     float inv = 1.0f / det;
//                     uval = (-S_Iy2 * S_IxIt + S_Ixy * S_IyIt) * inv;
//                     vval = ( S_Ixy * S_IxIt - S_Ix2 * S_IyIt) * inv;
//                 }
//                 unsigned char fb[4];
//                 write_float_le_bytes(fb, uval);
//                 send_bytes_bare(fb,4);
//                 /* store v later in separate loop? Simpler: stream u for whole frame first, then v.
//                    But to avoid large RAM, we stream u immediately; we still need to stream v too -
//                    we can't compute v separately without storing; but we can stream u and v interleaved.
//                    The code currently sends u then immediately v. That's acceptable (consumer must match).
//                 */
//                 write_float_le_bytes(fb, vval);
//                 send_bytes_bare(fb,4);
//             }
//         }
//         /* After finishing frame f, continue to next frame. */
//     }

//     /* finished */
//     const char *done = "LK bare-metal done\n";
//     for(const char *s = done; *s; ++s) sbicon_putc(*s);
//     while(1) asm volatile("wfi");
//     return 0;
// }
// #endif
