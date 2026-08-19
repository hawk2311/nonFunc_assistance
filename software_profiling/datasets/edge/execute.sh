#!/bin/bash

echo Starting 
python3 image.py
sleep 1
riscv64-unknown-elf-gcc -static  -o edgeBare.elf edgeBare.c  -lm
sleep 1
spike pk edgeBare.elf
