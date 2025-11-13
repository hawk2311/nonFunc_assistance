#ifndef CLOCKTIMER_H
#define CLOCKTIMER_H

/* Bare-metal version of ClockTimer.
 *
 * This implementation does NOT use gettimeofday (no syscalls). Instead it
 * converts cycle counts into microseconds using CPU_FREQ_HZ which must be
 * provided at compile time (as a macro or -D CPU_FREQ_HZ=...).
 *
 * Usage: compile with -DCPU_FREQ_HZ=50000000  (for 50 MHz) or define it
 * in a platform header.
 *
 * If you prefer to use a different timebase, replace getticks()/conversion.
 */

#include "cycle.h"
#include <stdint.h>

#ifndef CPU_FREQ_HZ
#warning "CPU_FREQ_HZ not defined — ClockTimer will not compile without CPU_FREQ_HZ. Define -DCPU_FREQ_HZ=<hz>."
#endif

#define CPU_FREQ_HZ 50000000ULL  //Need to be changed to the Rocket Chip Freq

class ClockTimer
{
public:
    void start()
    {
        t0 = getticks();
    }

    /* stop() returns microseconds as signed long (64-bit safe if using long long).
       For portability use long long return type if you want larger ranges. */
    long long stop()
    {
        ticks t1 = getticks();
        ticks diff = (t1 >= t0) ? (t1 - t0) : 0;
        /* microseconds = cycles / (CPU_FREQ_HZ / 1e6) = (cycles * 1_000_000) / CPU_FREQ_HZ */
        /* do multiplication first in 128-bit if available, otherwise 64-bit may overflow for very long intervals.
           We'll use unsigned long long and assume intervals are not huge in typical embedded use. */
#if defined(__SIZEOF_INT128__)
        __uint128_t tmp = ( __uint128_t ) diff * 1000000ULL;
        tmp /= ( __uint128_t ) CPU_FREQ_HZ;
        return (long long) tmp;
#else
        /* Fallback: do 64-bit math carefully. This will overflow if diff*1e6 > 2^64-1.
           For most short measurements this is OK. If you need very long intervals, enable __int128 or split math. */
        unsigned long long tmp = diff;
        tmp = tmp / (CPU_FREQ_HZ / 1000000ULL); /* integer division — loses precision */
        return (long long) tmp;
#endif
    }

private:
    ticks t0;
};

#endif /* CLOCKTIMER_H */

