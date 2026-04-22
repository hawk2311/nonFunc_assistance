#!/bin/bash

riscv64-unknown-elf-gcc -static -o optflow_2frames.elf optflow_2frames.cpp -lm

spike pk optflow_2frames.elf > flow_output.txt
echo "Spike simulation completed"
python3 create_image.py
