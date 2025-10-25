/*
 
 Einfacher, selbständiger Horn-Schunck dense optical flow
 Eingabe: zwei PGM (P5) Graustufenbilder gleicher Größe
 Ausgabe: eine PPM (P6) Farbbild visualisiert den Flow (HSV->RGB)
 Kein OpenCV, nur stdio / stdlib / math / string.
 */

 //Source: https://chatgpt.com/share/68fd0aa3-66d8-8005-a371-56e2f3278bd6
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef unsigned char u8;

static void die(const char *msg){
    fprintf(stderr, "%s\n", msg);
    exit(1);
}

/* --- PGM P5 lesen --- */
u8 *read_pgm(const char *fname, int *w, int *h){
    FILE *f = fopen(fname, "rb");
    if(!f) { perror("fopen"); die("Fehler beim Öffnen der PGM-Datei"); }
    char magic[3] = {0};
    if(fscanf(f, "%2s", magic)!=1) die("Ungültiges PGM");
    if(strcmp(magic, "P5")!=0) die("Nur P5 (binary) PGM unterstützt");
    // skip whitespace and comments
    int c;
    c = fgetc(f);
    while(c=='#'){ // comment line
        while(c!='\n') c = fgetc(f);
        c = fgetc(f);
    }
    ungetc(c, f);
    int width, height, maxv;
    if(fscanf(f, "%d %d %d", &width, &height, &maxv)!=3) die("Fehler beim Header lesen");
    if(maxv != 255) die("Nur 8-bit PGM unterstützt (max 255)");
    fgetc(f); // single whitespace before binary
    u8 *buf = malloc(width * height);
    if(!buf) die("Out of memory");
    size_t read = fread(buf, 1, width*height, f);
    if(read != (size_t)width*height) die("Fehler beim Lesen der Bilddaten");
    fclose(f);
    *w = width; *h = height;
    return buf;
}

/* --- PPM P6 schreiben --- */
void write_ppm(const char *fname, const unsigned char *rgb, int w, int h){
    FILE *f = fopen(fname, "wb");
    if(!f) { perror("fopen"); die("Fehler beim Öffnen der PPM-Datei"); }
    fprintf(f, "P6\n%d %d\n255\n", w, h);
    size_t written = fwrite(rgb, 1, (size_t)w*h*3, f);
    if(written != (size_t)w*h*3) die("Fehler beim Schreiben der PPM-Daten");
    fclose(f);
}

/* --- HSV (h in [0,360), s,v in [0,1]) -> RGB [0..255] --- */
void hsv_to_rgb(float h, float s, float v, unsigned char *r, unsigned char *g, unsigned char *b){
    float c = v * s;
    float hh = h / 60.0f;
    float x = c * (1 - fabsf(fmodf(hh, 2.0f) - 1));
    float m = v - c;
    float rr=0, gg=0, bb=0;
    if(hh >= 0 && hh < 1){ rr=c; gg=x; bb=0; }
    else if(hh < 2){ rr=x; gg=c; bb=0; }
    else if(hh < 3){ rr=0; gg=c; bb=x; }
    else if(hh < 4){ rr=0; gg=x; bb=c; }
    else if(hh < 5){ rr=x; gg=0; bb=c; }
    else { rr=c; gg=0; bb=x; }
    *r = (unsigned char)fminf(255.0f, (rr + m) * 255.0f);
    *g = (unsigned char)fminf(255.0f, (gg + m) * 255.0f);
    *b = (unsigned char)fminf(255.0f, (bb + m) * 255.0f);
}

/* --- Horn-Schunck dense optical flow --- 
   I1, I2: float images (w*h), output u,v (floats)
   alpha: smoothness weight
   num_iter: iterations
*/
void horn_schunck(const float *I1, const float *I2, float *u, float *v, int w, int h, float alpha, int num_iter){
    int N = w*h;
    float *Ix = malloc(sizeof(float)*N);
    float *Iy = malloc(sizeof(float)*N);
    float *It = malloc(sizeof(float)*N);
    if(!Ix || !Iy || !It) die("malloc failed");

    // compute gradients with simple finite differences
    for(int y=0;y<h;y++){
        for(int x=0;x<w;x++){
            int idx = y*w + x;
            int xm = x==0 ? x : x-1;
            int xp = x==w-1 ? x : x+1;
            int ym = y==0 ? y : y-1;
            int yp = y==h-1 ? y : y+1;
            float dx = (I1[y*w + xp] - I1[y*w + xm] + I2[y*w + xp] - I2[y*w + xm]) * 0.25f;
            float dy = (I1[yp*w + x] - I1[ym*w + x] + I2[yp*w + x] - I2[ym*w + x]) * 0.25f;
            float dt = (I2[idx] - I1[idx]);
            Ix[idx] = dx;
            Iy[idx] = dy;
            It[idx] = dt;
            u[idx] = 0.0f;
            v[idx] = 0.0f;
        }
    }

    // iterative update
    float alpha2 = alpha * alpha;
    float *u_avg = malloc(sizeof(float)*N);
    float *v_avg = malloc(sizeof(float)*N);
    if(!u_avg || !v_avg) die("malloc failed 2");

    for(int it=0; it<num_iter; ++it){
        // compute neighbor average (4-neighborhood)
        for(int y=0;y<h;y++){
            for(int x=0;x<w;x++){
                int idx = y*w + x;
                float sumu = 0.0f, sumv = 0.0f;
                int count = 0;
                if(x>0){ sumu += u[idx-1]; sumv += v[idx-1]; count++; }
                if(x<w-1){ sumu += u[idx+1]; sumv += v[idx+1]; count++; }
                if(y>0){ sumu += u[idx-w]; sumv += v[idx-w]; count++; }
                if(y<h-1){ sumu += u[idx+w]; sumv += v[idx+w]; count++; }
                u_avg[idx] = sumu / (count?count:1);
                v_avg[idx] = sumv / (count?count:1);
            }
        }
        // update u,v
        for(int i=0;i<N;i++){
            float ix = Ix[i], iy = Iy[i], itf = It[i];
            float num = ix * u_avg[i] + iy * v_avg[i] + itf;
            float den = alpha2 + ix*ix + iy*iy;
            float tmp = (den==0.0f) ? 0.0f : (num / den);
            u[i] = u_avg[i] - ix * tmp;
            v[i] = v_avg[i] - iy * tmp;
        }
    }

    free(Ix); free(Iy); free(It); free(u_avg); free(v_avg);
}

/* --- convert u,v flow to color image ---
   We map angle -> hue (0..360), magnitude -> value (0..1 with clipping)
*/
unsigned char *flow_to_color(const float *u, const float *v, int w, int h, float max_magnitude){
    int N = w*h;
    unsigned char *rgb = malloc(3*N);
    if(!rgb) die("malloc failed flow_to_color");
    for(int i=0;i<N;i++){
        float ux = u[i], vy = v[i];
        float ang = atan2f(vy, ux); // -pi .. pi
        float mag = sqrtf(ux*ux + vy*vy);
        float hue = (ang + M_PI) * (180.0f / M_PI); // 0..360
        float val = mag / (max_magnitude>0 ? max_magnitude : 1.0f);
        if(val>1.0f) val = 1.0f;
        float sat = 1.0f;
        unsigned char r,g,b;
        hsv_to_rgb(hue, sat, val, &r, &g, &b);
        rgb[3*i+0]=r; rgb[3*i+1]=g; rgb[3*i+2]=b;
    }
    return rgb;
}

/* --- helper: convert u8 image to float normalized [0..1] --- */
void u8_to_float(const u8 *in, float *out, int N){
    for(int i=0;i<N;i++) out[i] = (float)in[i] / 255.0f;
}

/* --- main: usage:
   horn_schunk <frame1.pgm> <frame2.pgm> <out.ppm> [alpha] [iter] [max_magnitude]
*/
int main(int argc, char **argv){
    if(argc < 4){
        fprintf(stderr, "Usage: %s frame1.pgm frame2.pgm out.ppm [alpha (default 1.0)] [iter (default 100)] [max_magnitude (default 1.0)]\n", argv[0]);
        return 1;
    }
    const char *f1 = argv[1];
    const char *f2 = argv[2];
    const char *out = argv[3];
    float alpha = (argc>4) ? atof(argv[4]) : 1.0f;
    int iter = (argc>5) ? atoi(argv[5]) : 100;
    float max_m = (argc>6) ? atof(argv[6]) : 1.0f;

    int w1,h1,w2,h2;
    u8 *b1 = read_pgm(f1, &w1, &h1);
    u8 *b2 = read_pgm(f2, &w2, &h2);
    if(w1!=w2 || h1!=h2) die("Bildgrößen stimmen nicht überein");
    int w = w1, h = h1;
    int N = w*h;

    float *I1 = malloc(sizeof(float)*N);
    float *I2 = malloc(sizeof(float)*N);
    float *u = malloc(sizeof(float)*N);
    float *v = malloc(sizeof(float)*N);
    if(!I1 || !I2 || !u || !v) die("malloc failed main");

    u8_to_float(b1, I1, N);
    u8_to_float(b2, I2, N);

    horn_schunck(I1, I2, u, v, w, h, alpha, iter);

    // compute a reasonable max magnitude if user left default 1.0
    if(max_m <= 0.0f){
        float mm = 0.0f;
        for(int i=0;i<N;i++){ float m = sqrtf(u[i]*u[i] + v[i]*v[i]); if(m>mm) mm = m; }
        max_m = mm>1e-6f ? mm : 1.0f;
    }

    unsigned char *rgb = flow_to_color(u, v, w, h, max_m);
    write_ppm(out, rgb, w, h);

    free(b1); free(b2); free(I1); free(I2); free(u); free(v); free(rgb);
    return 0;
}
