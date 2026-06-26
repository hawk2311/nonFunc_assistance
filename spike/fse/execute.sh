#!/bin/bash

echo "The Conda envirnoment must be activated, therefore execute the env.sh file in the chipyard directory!"

riscv64-unknown-elf-gcc -static -o fse.elf fse.c  -lm

spike pk fse.elf > output.txt 

./generate_image

echo "Done, reconstructed image generated in this directory."

