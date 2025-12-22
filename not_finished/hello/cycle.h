#ifndef CYCLE_H
#define CYCLE_H

/* Minimal, bare-metal capable cycle counter for RISC-V.
 *
 * Provides:
 *   typedef unsigned long long ticks;
 *   ticks getticks(void);
 *   double elapsed(ticks t1, ticks t0);
 *
 * For rv64: read CSRR cycle (64bit).
 * For rv32: read cycleh/cycle loop to form 64-bit value.
 *
 * Requires GCC/Clang (inline asm).
 */

#include <stdint.h>

#if defined(__riscv)

typedef unsigned long long ticks;

#if defined(__riscv_xlen) && (__riscv_xlen == 64)

/* RV64: cycle CSR is 64-bit */
static inline ticks getticks(void)
{
    ticks val;
    asm volatile ("csrr %0, cycle" : "=r"(val));
    return val;
}

#elif defined(__riscv_xlen) && (__riscv_xlen == 32)

/* RV32: cycle is 64-bit split into cycle (low) and cycleh (high).
   Need to read cycleh, cycle, cycleh and loop until high didn't change. */
static inline ticks getticks(void)
{
    uint32_t hi0, lo, hi1;
    ticks ret;
    do {
        asm volatile ("csrr %0, cycleh" : "=r"(hi0));
        asm volatile ("csrr %0, cycle"  : "=r"(lo));
        asm volatile ("csrr %0, cycleh" : "=r"(hi1));
    } while (hi0 != hi1);
    ret = ((unsigned long long)hi1 << 32) | (unsigned long long)lo;
    return ret;
}

#else
/* Fallback: if __riscv_xlen is not defined, attempt a runtime check.
   Assume rv64 if sizeof(void*)==8. */
#if UINTPTR_MAX == 0xffffffffffffffffULL
static inline ticks getticks(void)
{
    ticks val;
    asm volatile ("csrr %0, cycle" : "=r"(val));
    return val;
}
#else
static inline ticks getticks(void)
{
    uint32_t hi0, lo, hi1;
    ticks ret;
    do {
        asm volatile ("csrr %0, cycleh" : "=r"(hi0));
        asm volatile ("csrr %0, cycle"  : "=r"(lo));
        asm volatile ("csrr %0, cycleh" : "=r"(hi1));
    } while (hi0 != hi1);
    ret = ((unsigned long long)hi1 << 32) | (unsigned long long)lo;
    return ret;
}
#endif
#endif /* __riscv_xlen */

static inline double elapsed(ticks t1, ticks t0)
{
    return (double)t1 - (double)t0;
}

#else
#error "cycle.h: Not a RISC-V target (define __riscv) or add support for your arch"
#endif /* __riscv */

#endif /* CYCLE_H */

