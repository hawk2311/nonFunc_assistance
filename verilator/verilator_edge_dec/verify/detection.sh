#!/bin/bash

echo Starting 
python3 image.py
sleep 8
riscv64-unknown-elf-gcc -static  -o edgeBare_for_spike.elf edgeBare_for_spike.c  -lm
sleep 5
spike pk edgeBare_for_spike.elf
