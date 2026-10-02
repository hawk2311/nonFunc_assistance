The create_dataset.py works the same as in the parent directory. Cachegrind always gives the same values, that's why there is only one run for every image.

Start script for edge detection:
python3 create_dataset.py "gcc -static -o edgeBare edgeBare.c -lm" edgeBare.c edgeBare edge_dec_data.csv 

Start script for fftw:
python3 create_dataset.py "gcc -static -o fftw fftw.c -lfftw3 -lm" fftw.c fftw fftw_data.csv 


