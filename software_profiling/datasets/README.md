Using script in general:
python3 create_dataset.py "whole_command_for_compiling" c_code name_of_compiled_code data_from_all_runs avg_data {-header}

When executing the script you need to add:
1. The command you want to compile your code with (with all flags etc)
2. the name of your code you want to compile
3. the name of the compiled code (the executable)
4. the name of the CSV file which should store the data of all runs
5. the name of the CSV file for the calculated averages
(optional) 6. add the "-header" flag to generate for all images in the "image" directory the corresponding header files


Start script for edge detection:
python3 create_dataset.py "gcc -static -o edgeBare edgeBare.c -lm" edgeBare.c edgeBare edge_dec_data.csv edge_dec_avg.csv 

Start script for fftw:
python3 create_dataset.py "gcc -static -o fftw fftw.c -lfftw3 -lm" fftw.c fftw fftw_data.csv fftw_avg.csv 

