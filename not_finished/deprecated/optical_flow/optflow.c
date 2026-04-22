#include <stdint.h>
#include <stdio.h>
//#include "optical_flow.h"
#include "video_data.h"

/*
 * Erwartet aus video_data.h:
 *
 * #define VIDEO_T ...
 * #define VIDEO_H ...
 * #define VIDEO_W ...
 * extern const uint8_t video_data[VIDEO_T][VIDEO_H][VIDEO_W];
 */



 typedef struct {
    int u;   // Bewegung in x-Richtung
    int v;   // Bewegung in y-Richtung
} Flow;


static inline int Ix(
    const uint8_t frame[VIDEO_H][VIDEO_W],
    int x, int y)
{
    return (int)frame[y][x + 1] - (int)frame[y][x - 1];
}

static inline int Iy(
    const uint8_t frame[VIDEO_H][VIDEO_W],
    int x, int y)
{
    return (int)frame[y + 1][x] - (int)frame[y - 1][x];
}

static inline int It(
    const uint8_t f1[VIDEO_H][VIDEO_W],
    const uint8_t f2[VIDEO_H][VIDEO_W],
    int x, int y)
{
    return (int)f2[y][x] - (int)f1[y][x];
}

Flow optical_flow_data[VIDEO_T - 1][VIDEO_H][VIDEO_W];


static Flow optical_flow_pixel(
    const uint8_t f1[VIDEO_H][VIDEO_W],
    const uint8_t f2[VIDEO_H][VIDEO_W],
    int x, int y)
{
    int sum_ix2 = 0;
    int sum_iy2 = 0;
    int sum_ixy = 0;
    int sum_ixt = 0;
    int sum_iyt = 0;

    for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {

            int ix = Ix(f1, x + dx, y + dy);
            int iy = Iy(f1, x + dx, y + dy);
            int it = It(f1, f2, x + dx, y + dy);

            sum_ix2 += ix * ix;
            sum_iy2 += iy * iy;
            sum_ixy += ix * iy;
            sum_ixt += ix * it;
            sum_iyt += iy * it;
        }
    }

    Flow flow = {0, 0};

    int det = sum_ix2 * sum_iy2 - sum_ixy * sum_ixy;

    if (det != 0) {
        flow.u = (-sum_iy2 * sum_ixt + sum_ixy * sum_iyt) / det;
        flow.v = ( sum_ixy * sum_ixt - sum_ix2 * sum_iyt) / det;
    }

    return flow;
}


int main(){
    //compute_optical_flow_video();
     for (int t = 0; t < VIDEO_T - 1; t++) {

        const uint8_t (*f1)[VIDEO_W] = video_data[t];
        const uint8_t (*f2)[VIDEO_W] = video_data[t + 1];

        for (int y = 1; y < VIDEO_H - 1; y++) {
            for (int x = 1; x < VIDEO_W - 1; x++) {

                Flow flow = optical_flow_pixel(f1, f2, x, y);
                //optical_flow_data[t][y][x] = flow;
                //printf("Frame %d, Pixel (%d, %d): Flow (u=%d, v=%d)\n", t, x, y, flow.u, flow.v);
                printf("%d %d %d %d %d\n", t, x, y, flow.u, flow.v);


                // if ((flow.u != 0) || (flow.v != 0)) {
                //     // Beispiel: Bewegung erkannt
                //     // (z. B. Zähler erhöhen, Flag setzen, etc.)
                // }
            }
        }


    }
    return 0;
}



// #include <stdint.h>
// #include <stdio.h>
// #include "video_data.h"
// #include "optical_flow.h"

// /*
//  * Erwartet aus video_data.h:
//  *
//  * #define VIDEO_T ...
//  * #define VIDEO_H ...
//  * #define VIDEO_W ...
//  * extern const uint8_t video_data[VIDEO_T][VIDEO_H][VIDEO_W];
//  */

// typedef struct {
//     int16_t u;   // Bewegung in x-Richtung
//     int16_t v;   // Bewegung in y-Richtung
// } Flow;

// /*
//  * Optical-Flow-Ergebnis:
//  * [Frame][y][x]
//  * Frame bezieht sich auf (t -> t+1)
//  */
// //tatic Flow optical_flow_data[VIDEO_T - 1][VIDEO_H][VIDEO_W];

// /* -------------------------------------------------- */
// /* Ableitungen                                        */
// /* -------------------------------------------------- */

// static inline int Ix(
//     const uint8_t frame[VIDEO_H][VIDEO_W],
//     int x, int y)
// {
//     return (int)frame[y][x + 1] - (int)frame[y][x - 1];
// }

// static inline int Iy(
//     const uint8_t frame[VIDEO_H][VIDEO_W],
//     int x, int y)
// {
//     return (int)frame[y + 1][x] - (int)frame[y - 1][x];
// }

// static inline int It(
//     const uint8_t f1[VIDEO_H][VIDEO_W],
//     const uint8_t f2[VIDEO_H][VIDEO_W],
//     int x, int y)
// {
//     return (int)f2[y][x] - (int)f1[y][x];
// }

// /* -------------------------------------------------- */
// /* Optical Flow für ein Pixel (Lucas–Kanade, 3×3)     */
// /* -------------------------------------------------- */

// static Flow optical_flow_pixel(
//     const uint8_t f1[VIDEO_H][VIDEO_W],
//     const uint8_t f2[VIDEO_H][VIDEO_W],
//     int x, int y)
// {
//     int sum_ix2 = 0;
//     int sum_iy2 = 0;
//     int sum_ixy = 0;
//     int sum_ixt = 0;
//     int sum_iyt = 0;

//     for (int dy = -1; dy <= 1; dy++) {
//         for (int dx = -1; dx <= 1; dx++) {

//             int ix = Ix(f1, x + dx, y + dy);
//             int iy = Iy(f1, x + dx, y + dy);
//             int it = It(f1, f2, x + dx, y + dy);

//             sum_ix2 += ix * ix;
//             sum_iy2 += iy * iy;
//             sum_ixy += ix * iy;
//             sum_ixt += ix * it;
//             sum_iyt += iy * it;
//         }
//     }

//     Flow flow = {0, 0};

//     int det = sum_ix2 * sum_iy2 - sum_ixy * sum_ixy;

//     if (det != 0) {
//         flow.u = (int16_t)((-sum_iy2 * sum_ixt + sum_ixy * sum_iyt) / det);
//         flow.v = (int16_t)(( sum_ixy * sum_ixt - sum_ix2 * sum_iyt) / det);
//     }

//     return flow;
// }

// /* -------------------------------------------------- */
// /* Main                                               */
// /* -------------------------------------------------- */

// int main(void)
// {
//     /* Initialisierung (Sicherheitsmaßnahme) */
//     for (int t = 0; t < VIDEO_T - 1; t++) {
//         for (int y = 0; y < VIDEO_H; y++) {
//             for (int x = 0; x < VIDEO_W; x++) {
//                 optical_flow_data[t][y][x].u = 0;
//                 optical_flow_data[t][y][x].v = 0;
//             }
//         }
//     }

//     /* Optical Flow berechnen */
//     for (int t = 0; t < VIDEO_T - 1; t++) {

//         const uint8_t (*f1)[VIDEO_W] = video_data[t];
//         const uint8_t (*f2)[VIDEO_W] = video_data[t + 1];

//         for (int y = 1; y < VIDEO_H - 1; y++) {
//             for (int x = 1; x < VIDEO_W - 1; x++) {

//                 optical_flow_data[t][y][x] =
//                     optical_flow_pixel(f1, f2, x, y);
//             }
//         }
//     }

//     /*
//      * Ab hier:
//      * optical_flow_data[t][y][x].u / .v
//      * → extern auslesen & visualisieren
//      */

//     return 0;
// }
