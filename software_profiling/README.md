

perf:

commands:
- compile C code: gcc -o edge edgeBare.c -lm
- compile with riscv: riscv64-unknown-elf-gcc edgeBare.c -o edgeBare -lm
- create record: perf record ./edge (-g flag can bea added)
- call report: perf report

- using perf stat: perf stat -e instructions,cycles,branches,bus-cycles,ref-cycles,cache-references ./edge

change settings for access (error occurs without change):
sudo sysctl kernel.perf_event_paranoid=-1

notes:
- egde was compiled with the gcc compiler whereas edgeBare was compiled with riscv64-unknown-elf-gcc but it seems that there is an error when this is used in perf
- 


valgrind:

commands:
- using valgrind: valgrind --tool=callgrind ./edge
- kcachegrind callgrind.out.<pid>

