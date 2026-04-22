This implmentation calculates the  Optical Flow for the whole video and calculates the average movement in the whole scene. This average value is drawed in an output image, which is the last frame of the video with a line showing the average movement of the scene.

Workflow:

-Compiling the code for RISC V
riscv64-unknown-elf-gcc -static -o optflow.elf optflow.cpp -lm

-Execute the binary with spike and move the output to a txt file
spike pk optflow.elf > flow_output.txt

-with the help of the python script it is possible to generate an output image based on the values from the txt file
python3 create_image.py
