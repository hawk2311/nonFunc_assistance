# 0 "start.S"
# 0 "<built-in>"
# 0 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 0 "<command-line>" 2
# 1 "start.S"

.section .text.start
.global _start
_start:
    la sp, _stack_top
    call main
1: j 1b
