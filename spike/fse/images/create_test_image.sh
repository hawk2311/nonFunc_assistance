#!/bin/bash

echo "Creating test image for Frequency Selective Extrapolation."

#Generating an array in a C Header file which contains image data.
echo "Enter an image for your test purpose:"
python3 image.py
echo "Generated image data in header file."

#Generating an array in a C Header file which contains data of a selected mask which represent the errors in the image.
echo "Enter an mask (also as a image) which contains errors for your image purpose:"
python3 image.py
echo "Generated mask data in header file."

#This code will generate a overlay of the image and the mask and saves it in a new header file.
./create_masked_image

sleep 1

#This will show the overlay, this is just for validation.
./show_masked_image

#Conversion from pgm to jpg if needed. Probably useful for futher use in other code.
convert output.pgm output.jpg
