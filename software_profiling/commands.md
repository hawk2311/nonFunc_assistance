perf:

- compile C code: gcc -g -o edge edgeBare.c -lm
- create record: perf record ./edge (-g flag can bea added)
- call report: perf report
