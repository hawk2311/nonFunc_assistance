


//#include <stdio.h>
#include <stdint.h>
//#include <string.h>
#include <math.h>
#include "image_data.h"

#define WIDTH  133
#define HEIGHT 100
#define NPIX   (WIDTH * HEIGHT)

#define UART0_BASE 0x10000000

#define REG(base, offset) ((*((volatile unsigned char *)(base + offset))))
#define UART0_DR    REG(UART0_BASE, 0x00)
#define UART0_FCR   REG(UART0_BASE, 0x02)
#define UART0_LSR   REG(UART0_BASE, 0x05)
																						
#define UARTFCR_FFENA 0x01                // UART FIFO Control Register enable bit
#define UARTLSR_THRE 0x20                 // UART Line Status Register Transmit Hold Register Empty bit
#define UART0_FF_THR_EMPTY (UART0_LSR & UARTLSR_THRE)

static uint8_t gray[NPIX];
static uint8_t edges[NPIX];
static uint8_t overlay[NPIX * 3];   // RGB Overlay




void uart_putc(char c) {
  while (!UART0_FF_THR_EMPTY);            // Wait until the FIFO holding register is empty
  UART0_DR = c;                           // Write character to transmitter register
}

void uart_puts(char *str) {
  while (*str) {                          // Loop until value at string pointer is zero
    uart_putc(*str++);                    // Write the character and increment pointer
  }
}

int main() {

    //printf("[INFO] Starting edge detection...\n");
    uart_puts("Start");
    uart_puts("\n");
    uart_puts("\r");

    //----------------------------------------------------------------------
    // 1. Bild in eindimensionales Array kopieren
    //----------------------------------------------------------------------
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            //gray[y * WIDTH + x] = image_data[y][x];
            gray[y * WIDTH + x] = image_data[y][x];
        }
    }

    
    //----------------------------------------------------------------------
    // 2. Sobel-Kernel definieren
    //----------------------------------------------------------------------
    const int Gx[3][3] = {
        {-1, 0, 1},
        {-2, 0, 2},
        {-1, 0, 1}
    };
    const int Gy[3][3] = {
        {-1, -2, -1},
         {0,  0,  0},
         {1,  2,  1}
    };

    //----------------------------------------------------------------------
    // 3. Sobel-Kantenberechnung
    //----------------------------------------------------------------------
    for (int y = 1; y < HEIGHT - 1; y++) {
        for (int x = 1; x < WIDTH - 1; x++) {

            int sumX = 0;
            int sumY = 0;

            for (int ky = -1; ky <= 1; ky++) {
                for (int kx = -1; kx <= 1; kx++) {
                    int pixel = gray[(y + ky) * WIDTH + (x + kx)];
                    sumX += pixel * Gx[ky + 1][kx + 1];
                    sumY += pixel * Gy[ky + 1][kx + 1];
                }
            }
            
            //int magnitude = (int)sqrt((double)(sumX * sumX + sumY * sumY));
            char magnitude = 0;
            if (magnitude > 255) magnitude = 255;
            
            //edges[y * WIDTH + x] = (uint8_t)magnitude;
            uart_puts("Checkpoint 6");
            uart_puts("\n");
            uart_puts("\r");

            int new_value[4]; 	
            char buf_temp[4]; 	
    	 	
            for (int x = 0; x < 4; x++){ 		
                new_value[x] = magnitude % 10; 		
                magnitude = magnitude / 10; 		
                buf_temp[x] = new_value[x] +48;
            }

            uart_puts();
        }
    }

    //----------------------------------------------------------------------
    // 4. Overlay-Bild erzeugen (reines RGB-Array)
    //----------------------------------------------------------------------
    for (int i = 0; i < NPIX; i++) {
        uint8_t v = gray[i];

        overlay[3*i + 0] = v;
        overlay[3*i + 1] = v;
        overlay[3*i + 2] = v;

        if (edges[i] > 100) {
            overlay[3*i + 0] = 255; // rot
            overlay[3*i + 1] = 0;
            overlay[3*i + 2] = 0;
        }
    }

    //----------------------------------------------------------------------
    // 5. Ergebnis in export-Array kopieren (für Host)
    //----------------------------------------------------------------------
    // for (int i = 0; i < NPIX; i++)
    //     edges_export[i] = edges[i];

    
    //for (int y = 0; y < HEIGHT; y++) {
    //    for (int x = 0; x < WIDTH; x++) {
    //        //("edges: %i\n", edges[y * WIDTH + x]);
    //        //edges_output[y][x] = edges[y * WIDTH + x];
    //        printf("%i\n", edges[y * WIDTH + x]);
    //    }
    //}

    //----------------------------------------------------------------------
    // 6. Kontrollsumme ausgeben
    //----------------------------------------------------------------------
    // unsigned long sum = 0;
    // for (int i = 0; i < NPIX; i++) sum += edges[i];

    //printf("[INFO] Average edge intensity = %lu\n", sum / NPIX);
    //printf("[INFO] Done.\n");
    uart_puts("Done");

    return 0;
}



