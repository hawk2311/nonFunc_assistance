#!/bin/bash

python3 create_dataset.py "gcc -O2 -static -o edgeBare edgeBare.c -lm" edgeBare.c edgeBare edge_dec_data.csv edge_dec_avg.csv 

python3 create_dataset.py "gcc -O2 -static -o fftw fftw.c -lfftw3 -lm" fftw.c fftw fftw_data.csv fftw_avg.csv 

python3 merge_csv.py "edge_dec_avg.csv","fftw_avg.csv" all_data.csv