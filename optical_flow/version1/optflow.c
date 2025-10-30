// /*
 
//  Einfacher, selbständiger Horn-Schunck dense optical flow
//  Eingabe: zwei PGM (P5) Graustufenbilder gleicher Größe
//  Ausgabe: eine PPM (P6) Farbbild visualisiert den Flow (HSV->RGB)
//  Kein OpenCV, nur stdio / stdlib / math / string.
//  */

//  //Source: https://chatgpt.com/share/68fd0aa3-66d8-8005-a371-56e2f3278bd6
 
// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>
// #include <string.h>

// #ifndef M_PI
// #define M_PI 3.14159265358979323846
// #endif

// typedef unsigned char u8;

// static void die(const char *msg){
//     fprintf(stderr, "%s\n", msg);
//     exit(1);
// }

// /* --- PGM P5 lesen --- */
// u8 *read_pgm(const char *fname, int *w, int *h){
//     FILE *f = fopen(fname, "rb");
//     if(!f) { perror("fopen"); die("Fehler beim Öffnen der PGM-Datei"); }
//     char magic[3] = {0};
//     if(fscanf(f, "%2s", magic)!=1) die("Ungültiges PGM");
//     if(strcmp(magic, "P5")!=0) die("Nur P5 (binary) PGM unterstützt");
//     // skip whitespace and comments
//     int c;
//     c = fgetc(f);
//     while(c=='#'){ // comment line
//         while(c!='\n') c = fgetc(f);
//         c = fgetc(f);
//     }
//     ungetc(c, f);
//     int width, height, maxv;
//     if(fscanf(f, "%d %d %d", &width, &height, &maxv)!=3) die("Fehler beim Header lesen");
//     if(maxv != 255) die("Nur 8-bit PGM unterstützt (max 255)");
//     fgetc(f); // single whitespace before binary
//     u8 *buf = malloc(width * height);
//     if(!buf) die("Out of memory");
//     size_t read = fread(buf, 1, width*height, f);
//     if(read != (size_t)width*height) die("Fehler beim Lesen der Bilddaten");
//     fclose(f);
//     *w = width; *h = height;
//     return buf;
// }

// /* --- PPM P6 schreiben --- */
// void write_ppm(const char *fname, const unsigned char *rgb, int w, int h){
//     FILE *f = fopen(fname, "wb");
//     if(!f) { perror("fopen"); die("Fehler beim Öffnen der PPM-Datei"); }
//     fprintf(f, "P6\n%d %d\n255\n", w, h);
//     size_t written = fwrite(rgb, 1, (size_t)w*h*3, f);
//     if(written != (size_t)w*h*3) die("Fehler beim Schreiben der PPM-Daten");
//     fclose(f);
// }

// /* --- HSV (h in [0,360), s,v in [0,1]) -> RGB [0..255] --- */
// void hsv_to_rgb(float h, float s, float v, unsigned char *r, unsigned char *g, unsigned char *b){
//     float c = v * s;
//     float hh = h / 60.0f;
//     float x = c * (1 - fabsf(fmodf(hh, 2.0f) - 1));
//     float m = v - c;
//     float rr=0, gg=0, bb=0;
//     if(hh >= 0 && hh < 1){ rr=c; gg=x; bb=0; }
//     else if(hh < 2){ rr=x; gg=c; bb=0; }
//     else if(hh < 3){ rr=0; gg=c; bb=x; }
//     else if(hh < 4){ rr=0; gg=x; bb=c; }
//     else if(hh < 5){ rr=x; gg=0; bb=c; }
//     else { rr=c; gg=0; bb=x; }
//     *r = (unsigned char)fminf(255.0f, (rr + m) * 255.0f);
//     *g = (unsigned char)fminf(255.0f, (gg + m) * 255.0f);
//     *b = (unsigned char)fminf(255.0f, (bb + m) * 255.0f);
// }

// /* --- Horn-Schunck dense optical flow --- 
//    I1, I2: float images (w*h), output u,v (floats)
//    alpha: smoothness weight
//    num_iter: iterations
// */
// void horn_schunck(const float *I1, const float *I2, float *u, float *v, int w, int h, float alpha, int num_iter){
//     int N = w*h;
//     float *Ix = malloc(sizeof(float)*N);
//     float *Iy = malloc(sizeof(float)*N);
//     float *It = malloc(sizeof(float)*N);
//     if(!Ix || !Iy || !It) die("malloc failed");

//     // compute gradients with simple finite differences
//     for(int y=0;y<h;y++){
//         for(int x=0;x<w;x++){
//             int idx = y*w + x;
//             int xm = x==0 ? x : x-1;
//             int xp = x==w-1 ? x : x+1;
//             int ym = y==0 ? y : y-1;
//             int yp = y==h-1 ? y : y+1;
//             float dx = (I1[y*w + xp] - I1[y*w + xm] + I2[y*w + xp] - I2[y*w + xm]) * 0.25f;
//             float dy = (I1[yp*w + x] - I1[ym*w + x] + I2[yp*w + x] - I2[ym*w + x]) * 0.25f;
//             float dt = (I2[idx] - I1[idx]);
//             Ix[idx] = dx;
//             Iy[idx] = dy;
//             It[idx] = dt;
//             u[idx] = 0.0f;
//             v[idx] = 0.0f;
//         }
//     }

//     // iterative update
//     float alpha2 = alpha * alpha;
//     float *u_avg = malloc(sizeof(float)*N);
//     float *v_avg = malloc(sizeof(float)*N);
//     if(!u_avg || !v_avg) die("malloc failed 2");

//     for(int it=0; it<num_iter; ++it){
//         // compute neighbor average (4-neighborhood)
//         for(int y=0;y<h;y++){
//             for(int x=0;x<w;x++){
//                 int idx = y*w + x;
//                 float sumu = 0.0f, sumv = 0.0f;
//                 int count = 0;
//                 if(x>0){ sumu += u[idx-1]; sumv += v[idx-1]; count++; }
//                 if(x<w-1){ sumu += u[idx+1]; sumv += v[idx+1]; count++; }
//                 if(y>0){ sumu += u[idx-w]; sumv += v[idx-w]; count++; }
//                 if(y<h-1){ sumu += u[idx+w]; sumv += v[idx+w]; count++; }
//                 u_avg[idx] = sumu / (count?count:1);
//                 v_avg[idx] = sumv / (count?count:1);
//             }
//         }
//         // update u,v
//         for(int i=0;i<N;i++){
//             float ix = Ix[i], iy = Iy[i], itf = It[i];
//             float num = ix * u_avg[i] + iy * v_avg[i] + itf;
//             float den = alpha2 + ix*ix + iy*iy;
//             float tmp = (den==0.0f) ? 0.0f : (num / den);
//             u[i] = u_avg[i] - ix * tmp;
//             v[i] = v_avg[i] - iy * tmp;
//         }
//     }

//     free(Ix); free(Iy); free(It); free(u_avg); free(v_avg);
// }

// /* --- convert u,v flow to color image ---
//    We map angle -> hue (0..360), magnitude -> value (0..1 with clipping)
// */
// unsigned char *flow_to_color(const float *u, const float *v, int w, int h, float max_magnitude){
//     int N = w*h;
//     unsigned char *rgb = malloc(3*N);
//     if(!rgb) die("malloc failed flow_to_color");
//     for(int i=0;i<N;i++){
//         float ux = u[i], vy = v[i];
//         float ang = atan2f(vy, ux); // -pi .. pi
//         float mag = sqrtf(ux*ux + vy*vy);
//         float hue = (ang + M_PI) * (180.0f / M_PI); // 0..360
//         float val = mag / (max_magnitude>0 ? max_magnitude : 1.0f);
//         if(val>1.0f) val = 1.0f;
//         float sat = 1.0f;
//         unsigned char r,g,b;
//         hsv_to_rgb(hue, sat, val, &r, &g, &b);
//         rgb[3*i+0]=r; rgb[3*i+1]=g; rgb[3*i+2]=b;
//     }
//     return rgb;
// }

// /* --- helper: convert u8 image to float normalized [0..1] --- */
// void u8_to_float(const u8 *in, float *out, int N){
//     for(int i=0;i<N;i++) out[i] = (float)in[i] / 255.0f;
// }

// /* --- main: usage:
//    horn_schunk <frame1.pgm> <frame2.pgm> <out.ppm> [alpha] [iter] [max_magnitude]
// */
// int main(int argc, char **argv){
//     if(argc < 4){
//         fprintf(stderr, "Usage: %s frame1.pgm frame2.pgm out.ppm [alpha (default 1.0)] [iter (default 100)] [max_magnitude (default 1.0)]\n", argv[0]);
//         return 1;
//     }
//     const char *f1 = argv[1];
//     const char *f2 = argv[2];
//     const char *out = argv[3];
//     float alpha = (argc>4) ? atof(argv[4]) : 1.0f;
//     int iter = (argc>5) ? atoi(argv[5]) : 100;
//     float max_m = (argc>6) ? atof(argv[6]) : 1.0f;

//     int w1,h1,w2,h2;
//     u8 *b1 = read_pgm(f1, &w1, &h1);
//     u8 *b2 = read_pgm(f2, &w2, &h2);
//     if(w1!=w2 || h1!=h2) die("Bildgrößen stimmen nicht überein");
//     int w = w1, h = h1;
//     int N = w*h;

//     float *I1 = malloc(sizeof(float)*N);
//     float *I2 = malloc(sizeof(float)*N);
//     float *u = malloc(sizeof(float)*N);
//     float *v = malloc(sizeof(float)*N);
//     if(!I1 || !I2 || !u || !v) die("malloc failed main");

//     u8_to_float(b1, I1, N);
//     u8_to_float(b2, I2, N);

//     horn_schunck(I1, I2, u, v, w, h, alpha, iter);

//     // compute a reasonable max magnitude if user left default 1.0
//     if(max_m <= 0.0f){
//         float mm = 0.0f;
//         for(int i=0;i<N;i++){ float m = sqrtf(u[i]*u[i] + v[i]*v[i]); if(m>mm) mm = m; }
//         max_m = mm>1e-6f ? mm : 1.0f;
//     }

//     unsigned char *rgb = flow_to_color(u, v, w, h, max_m);
//     write_ppm(out, rgb, w, h);

//     free(b1); free(b2); free(I1); free(I2); free(u); free(v); free(rgb);
//     return 0;
// }

/*
 * Führt Optical Flow (Horn-Schunck) auf einem Video aus:
 *  1. Zerlegt das Video mit ffmpeg in Graustufen-PGM-Frames
 *  2. Berechnet Optical Flow zwischen jedem aufeinanderfolgenden Frame-Paar
 *  3. Erstellt aus den Flow-Bildern wieder ein Video
 *  4. Löscht die temporären Frames
 *
 *  Keine OpenCV- oder sonstige Bibliotheken außer C-Standard + ffmpeg.
 * Source: https://chatgpt.com/c/69007566-3acc-8328-bab0-517abf7c2605
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef unsigned char u8;

static void die(const char *msg) {
    fprintf(stderr, "%s\n", msg);
    exit(1);
}

/* --- PGM (P5) lesen --- */
u8 *read_pgm(const char *fname, int *w, int *h) {
    FILE *f = fopen(fname, "rb");
    if(!f) { perror(fname); die("Fehler beim Öffnen der PGM-Datei"); }

    char magic[3] = {0};
    if(fscanf(f, "%2s", magic)!=1 || strcmp(magic, "P5")!=0)
        die("Nur P5 (binary) PGM unterstützt");

    int c;
    c = fgetc(f);
    while(c=='#') { while(c!='\n' && c!=EOF) c=fgetc(f); c=fgetc(f); }
    ungetc(c, f);

    int width, height, maxv;
    if(fscanf(f, "%d %d %d", &width, &height, &maxv)!=3) die("Fehler im Header");
    if(maxv!=255) die("Nur 8-Bit-PGM unterstützt");
    fgetc(f);
    u8 *buf = malloc(width*height);
    fread(buf,1,width*height,f);
    fclose(f);
    *w = width; *h = height;
    return buf;
}

/* --- PPM schreiben --- */
void write_ppm(const char *fname, const unsigned char *rgb, int w, int h) {
    FILE *f = fopen(fname, "wb");
    if(!f) { perror(fname); die("Fehler beim Öffnen der PPM-Datei"); }
    fprintf(f, "P6\n%d %d\n255\n", w, h);
    fwrite(rgb,1,w*h*3,f);
    fclose(f);
}

/* --- HSV → RGB --- */
void hsv_to_rgb(float h, float s, float v, unsigned char *r, unsigned char *g, unsigned char *b) {
    float c = v * s;
    float hh = h / 60.0f;
    float x = c * (1 - fabsf(fmodf(hh, 2.0f) - 1));
    float m = v - c;
    float rr=0, gg=0, bb=0;
    if(hh>=0 && hh<1){ rr=c; gg=x; bb=0; }
    else if(hh<2){ rr=x; gg=c; bb=0; }
    else if(hh<3){ rr=0; gg=c; bb=x; }
    else if(hh<4){ rr=0; gg=x; bb=c; }
    else if(hh<5){ rr=x; gg=0; bb=c; }
    else { rr=c; gg=0; bb=x; }
    *r=(unsigned char)((rr+m)*255.0f);
    *g=(unsigned char)((gg+m)*255.0f);
    *b=(unsigned char)((bb+m)*255.0f);
}

/* --- Horn-Schunck Algorithmus --- */
void horn_schunck(const float *I1, const float *I2, float *u, float *v, int w, int h, float alpha, int num_iter) {
    int N = w*h;
    float *Ix = malloc(sizeof(float)*N);
    float *Iy = malloc(sizeof(float)*N);
    float *It = malloc(sizeof(float)*N);
    float *u_avg = malloc(sizeof(float)*N);
    float *v_avg = malloc(sizeof(float)*N);

    for(int y=0;y<h;y++){
        for(int x=0;x<w;x++){
            int idx=y*w+x;
            int xm=x==0?x:x-1, xp=x==w-1?x:x+1;
            int ym=y==0?y:y-1, yp=y==h-1?y:y+1;
            float dx=(I1[y*w+xp]-I1[y*w+xm]+I2[y*w+xp]-I2[y*w+xm])*0.25f;
            float dy=(I1[yp*w+x]-I1[ym*w+x]+I2[yp*w+x]-I2[ym*w+x])*0.25f;
            Ix[idx]=dx; Iy[idx]=dy; It[idx]=I2[idx]-I1[idx];
            u[idx]=v[idx]=0.0f;
        }
    }

    float a2 = alpha*alpha;
    for(int iter=0;iter<num_iter;iter++){
        for(int y=0;y<h;y++){
            for(int x=0;x<w;x++){
                int idx=y*w+x;
                float sumu=0,sumv=0; int c=0;
                if(x>0){sumu+=u[idx-1];sumv+=v[idx-1];c++;}
                if(x<w-1){sumu+=u[idx+1];sumv+=v[idx+1];c++;}
                if(y>0){sumu+=u[idx-w];sumv+=v[idx-w];c++;}
                if(y<h-1){sumu+=u[idx+w];sumv+=v[idx+w];c++;}
                u_avg[idx]=sumu/(c?c:1);
                v_avg[idx]=sumv/(c?c:1);
            }
        }
        for(int i=0;i<N;i++){
            float ix=Ix[i], iy=Iy[i], it=It[i];
            float num=ix*u_avg[i]+iy*v_avg[i]+it;
            float den=a2+ix*ix+iy*iy;
            float q=(den==0)?0:num/den;
            u[i]=u_avg[i]-ix*q;
            v[i]=v_avg[i]-iy*q;
        }
    }
    free(Ix); free(Iy); free(It); free(u_avg); free(v_avg);
}

/* --- Flow zu Farbe --- */
unsigned char *flow_to_color(const float *u, const float *v, int w, int h, float max_m) {
    unsigned char *rgb = malloc(3*w*h);
    for(int i=0;i<w*h;i++){
        float mag = sqrtf(u[i]*u[i]+v[i]*v[i]);
        float ang = atan2f(v[i],u[i]);
        float hue = (ang+M_PI)*180.0f/M_PI;
        float val = fminf(1.0f, mag/(max_m>0?max_m:1.0f));
        unsigned char r,g,b;
        hsv_to_rgb(hue,1.0f,val,&r,&g,&b);
        rgb[3*i]=r; rgb[3*i+1]=g; rgb[3*i+2]=b;
    }
    return rgb;
}

void u8_to_float(const u8 *in, float *out, int N){
    for(int i=0;i<N;i++) out[i] = (float)in[i]/255.0f;
}

/* --- Hauptprogramm --- */
int main(int argc, char **argv) {
    if(argc < 2) {
        printf("Usage: %s input_video.mp4\n", argv[0]);
        return 1;
    }

    const char *video = argv[1];
    const char *frames_dir = "frames";
    const char *flows_dir = "flows";

    // 1. Frames extrahieren
    printf("[*] Zerlege Video in Frames...\n");
    char cmd[512];
    snprintf(cmd, sizeof(cmd),
        "rm -rf %s %s; mkdir %s %s; "
        "ffmpeg -hide_banner -loglevel error -i \"%s\" -pix_fmt gray -f image2 %s/frame%%04d.pgm",
        frames_dir, flows_dir, frames_dir, flows_dir, video, frames_dir);
    if(system(cmd)!=0) die("Fehler: ffmpeg konnte das Video nicht zerlegen.");

    // Zähle Frames
    int frame_count = 0;
    while(1){
        char fname[256];
        snprintf(fname, sizeof(fname), "%s/frame%04d.pgm", frames_dir, frame_count+1);
        FILE *f = fopen(fname, "rb");
        if(!f) break;
        fclose(f);
        frame_count++;
    }
    if(frame_count < 2) die("Zu wenige Frames im Video.");

    printf("[*] %d Frames gefunden.\n", frame_count);

    // 2. Optical Flow berechnen
    for(int i=1;i<frame_count;i++){
        char f1[256], f2[256], outname[256];
        snprintf(f1,sizeof(f1),"%s/frame%04d.pgm",frames_dir,i);
        snprintf(f2,sizeof(f2),"%s/frame%04d.pgm",frames_dir,i+1);
        snprintf(outname,sizeof(outname),"%s/flow%04d.ppm",flows_dir,i);

        int w1,h1,w2,h2;
        u8 *b1 = read_pgm(f1,&w1,&h1);
        u8 *b2 = read_pgm(f2,&w2,&h2);
        if(w1!=w2 || h1!=h2) die("Framesize mismatch!");

        int N=w1*h1;
        float *I1=malloc(sizeof(float)*N);
        float *I2=malloc(sizeof(float)*N);
        float *u=malloc(sizeof(float)*N);
        float *v=malloc(sizeof(float)*N);
        u8_to_float(b1,I1,N);
        u8_to_float(b2,I2,N);

        horn_schunck(I1,I2,u,v,w1,h1,1.0f,100);

        unsigned char *rgb = flow_to_color(u,v,w1,h1,1.0f);
        write_ppm(outname,rgb,w1,h1);
        printf("[*] Flow %d/%d gespeichert -> %s\n",i,frame_count-1,outname);

        free(b1); free(b2); free(I1); free(I2); free(u); free(v); free(rgb);
    }

    // 3. Flow-Bilder wieder zu Video
    printf("[*] Erstelle Output-Video aus Flow-Frames...\n");
    snprintf(cmd,sizeof(cmd),
        "ffmpeg -hide_banner -loglevel error -framerate 25 -i %s/flow%%04d.ppm "
        "-pix_fmt yuv420p -y output_flow.mp4",
        flows_dir);
    if(system(cmd)!=0) die("Fehler beim Erstellen des Ausgabevideos.");

    // 4. Aufräumen
    printf("[*] Lösche temporäre Frames...\n");
    snprintf(cmd,sizeof(cmd),"rm -rf %s", frames_dir);
    system(cmd);

    printf("[✓] Fertig! Output: output_flow.mp4\n");
    return 0;
}
