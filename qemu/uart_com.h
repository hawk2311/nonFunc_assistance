#define UART0_BASE 0x10000000

#define REG(base, offset) ((*((volatile unsigned char *)(base + offset))))
#define UART0_DR    REG(UART0_BASE, 0x00)
#define UART0_FCR   REG(UART0_BASE, 0x02)
#define UART0_LSR   REG(UART0_BASE, 0x05)
																						
#define UARTFCR_FFENA 0x01                // UART FIFO Control Register enable bit
#define UARTLSR_THRE 0x20                 // UART Line Status Register Transmit Hold Register Empty bit
#define UART0_FF_THR_EMPTY (UART0_LSR & UARTLSR_THRE)



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
