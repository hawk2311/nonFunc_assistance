#!/bin/bash
riscv64-unknown-elf-gcc -march=rv64imafd -mabi=lp64d  -mcmodel=medany -ffreestanding -c boot1.S -o boot1.o
riscv64-unknown-elf-gcc -march=rv64imafd  -mabi=lp64d -mcmodel=medany -ffreestanding -c edgeDetect_bm_Sync.c -o edgeDetect_bm_Sync.o
riscv64-unknown-elf-gcc -march=rv64imafd  -mabi=lp64d -mcmodel=medany -ffreestanding -nostdlib   -T link_uart.ld boot1.o edgeDetect_bm_Sync.o -lm -lgcc -o edgeDetect_bm_Sync.elf
riscv64-unknown-elf-objcopy -O binary edgeDetect_bm_Sync.elf edgeDetect_bm_Sync.bin
