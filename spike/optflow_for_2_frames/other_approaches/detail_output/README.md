With this approach every movement in between two frames is shown in the result image.
The workflow is:

-Compiling the code for RISC V
riscv64-unknown-elf-gcc -static -o optflow_2frames.elf optflow_2frames.cpp -lm

-Execute the binary with spike and move the output to a txt file
spike pk optflow_2frames.elf > flow_output.txt

-with the help of the python script it is possible to generate an output image based on the values from the txt file
python3 create_image.py

