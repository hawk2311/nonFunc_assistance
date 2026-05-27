Date of documentation: 27.05.2026

To run the code in edgeBare.c in qemu you need to install qemu-system-riscv64. 

With the help of the Makefile you can start execution, for this type "make run" in the command line, this will print the result. To store the result use "make run > output.txt". Now you can open the "output.txt" and delete the first lines of text, you just need the numbers. After saving this you can execute "./generate_image" and it will give you the file "result.pgm" which shows the result of the calculation.
The used image can be seen in /nonFunc_repo/verilator/verilator_edge_dec/ with the name "test.jpg".
With "make clean" you can delete the files created by "make run".
The code in "my_math.h" provides a simple way to calculate the square root for a given number. This was implemented due to problems with linking general math libaries in qemu. This could change in future, for now the result of my_sqrt (used in edgeBare.c, implemented in my_math.h) is sufficient for its purpose.  
