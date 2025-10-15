#!/bin/bash

echo Starting 
python3 image.py
sleep 8
riscv64-unknown-elf-gcc -static  -o edge.elf edge.c  -lm
sleep 5
spike pk edge.elf
