# default: edgeBare.elf

# edgeBare.o: edgeBare.c
# 	riscv64-unknown-elf-gcc -c -g -O0 -ffreestanding -march=rv32i -mabi=ilp32 -o edgeBare.o edgeBare.c -lm -lgcc

# start.o: start.s
# 	riscv64-unknown-elf-as -march=rv32i -mabi=ilp32 -o start.o start.s

# edgeBare.elf: edgeBare.o start.o baremetal.ld
# 	riscv64-unknown-elf-ld -T baremetal.ld -m elf32lriscv -o edgeBare.elf edgeBare.o start.o 

# run: edgeBare.elf
# 	@echo "Ctrl-A C for QEMU console, then quit to exit"
# 	qemu-system-riscv32 -nographic -serial mon:stdio -machine virt -bios edgeBare.elf

# clean:
# 	rm -f edgeBare.o start.o edgeBare.elf

default: optflow_2frames.elf

optflow_2frames.o: optflow_2frames.c
	riscv64-unknown-elf-gcc -c -g -O0 -ffreestanding -march=rv64imafdc -mabi=lp64d -mcmodel=medany -o optflow_2frames.o optflow_2frames.c

start.o: start.s
	riscv64-unknown-elf-as -march=rv64imafdc -mabi=lp64d -o start.o start.s

optflow_2frames.elf: optflow_2frames.o start.o baremetal.ld
	riscv64-unknown-elf-gcc -T baremetal.ld -march=rv64imafdc -mabi=lp64d -mcmodel=medany -nostartfiles -o optflow_2frames.elf optflow_2frames.o start.o -lm -lgcc

run: optflow_2frames.elf
#	@echo "Ctrl-A C for QEMU console, then quit to exit"
	qemu-system-riscv64 -nographic -serial mon:stdio -machine virt -bios optflow_2frames.elf

clean:
	rm -f optflow_2frames.o start.o optflow_2frames.elf
