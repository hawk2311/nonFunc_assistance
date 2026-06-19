

// main.c  


#include <stdint.h>

#define UART0_BASE 0x10000000

#define UART_THR  (*(volatile uint8_t *)(UART0_BASE + 0x0))  // Transmit Holding Register
#define UART_LSR  (*(volatile uint8_t *)(UART0_BASE + 0x5))  // Line Status Register
#define UART_LSR_THRE  0x20                                   // THR Empty bit

void uart_putchar(char c) {
    // Warten bis der UART bereit ist zu senden
    while ((UART_LSR & UART_LSR_THRE) == 0);
    UART_THR = c;
}

int main() {
    const char *s = "Hello RISC-V!\n";
    while (*s) uart_putchar(*s++);

    // Endlosschleife am Ende (statt return)
    while (1);
    return 0;
}
