


//#include <stdio.h>
#include <stdint.h>
//#include <string.h>
//#include <math.h>
#include "my_math.h"
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

// convert interger to string
void uart_putint(int n) {
    if (n == 0) {
        uart_putc('0');
        return;
    }

    char buf[12];
    int i = 0;

    
    if (n < 0) {
        uart_putc('-');
        n = -n;
    }

    // write values backwards in buffer
    while (n > 0) {
        buf[i++] = '0' + (n % 10); // get the least significant number with modulo
        n /= 10; //move the comma to next position, cut last number
    }

    // print buffer in other direction
    for (int j = i - 1; j >= 0; j--) {
        uart_putc(buf[j]);
    }
}


int main() {

    //printf("[INFO] Starting edge detection...\n");
    //uart_puts("Start");
    // uart_puts("\n");
    // uart_puts("\r");

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
            
            int magnitude = (int)my_sqrt((sumX * sumX + sumY * sumY)); 
            //int magnitude = (int)sqrt((double)(sumX * sumX + sumY * sumY)); //PROBLEM
            if (magnitude > 255) magnitude = 255;
            edges[y * WIDTH + x] = (uint8_t)magnitude;
            
            
            
        }
    }


    
    for (int y = 0; y < HEIGHT; y++) {
       for (int x = 0; x < WIDTH; x++) {
           uart_putint(edges[y * WIDTH + x]);
           uart_puts("\n");
       }
    }

  

    //uart_puts("Done");

    return 0;
}



