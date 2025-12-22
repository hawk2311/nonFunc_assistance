#!/bin/bash

echo Starting 
python3 image.py
sleep 5
riscv64-unknown-elf-gcc -static   -I./fftw-3.3.10/install/include   -L./fftw-3.3.10/install/lib -o fftw.elf fftw.c -lfftw3 -lm
sleep 5
spike pk fftw.elf
sleep 10
./gen_freq_pic

