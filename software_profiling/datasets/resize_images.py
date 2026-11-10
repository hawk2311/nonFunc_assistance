import subprocess
import os
import re
import csv
import sys
from PIL import Image
import numpy as np

#convert -resize 120% zelda.ppm test.png


def check_size(image):
    img= Image.open("images/"+image).convert("L")
    size = img.size #img.size returns a tuple
    if size[0]<700 and size[1]<700: #this is the category of "small" images, they all get scaled uo
        subprocess.run(["convert", "-resize" , "120%", "images/"+image, "images/120_"+image])
        subprocess.run(["convert", "-resize" , "140%", "images/"+image, "images/140_"+image])
        subprocess.run(["convert", "-resize" , "160%", "images/"+image, "images/160_"+image])
    elif (size[0]>=700 and size[0]<1000) or (size[1]>=700 and size[1]<1000) : #catergory of "middle-sized" images
        subprocess.run(["convert", "-resize" , "80%", "images/"+image, "images/80_"+image])
        subprocess.run(["convert", "-resize" , "120%", "images/"+image, "images/120_"+image])
        subprocess.run(["convert", "-resize" , "140%", "images/"+image, "images/140_"+image])
    else: #the rest are large images
        subprocess.run(["convert", "-resize" , "80%", "images/"+image, "images/80_"+image])
        subprocess.run(["convert", "-resize" , "60%", "images/"+image, "images/60_"+image])
        subprocess.run(["convert", "-resize" , "40%", "images/"+image, "images/40_"+image])
    


def main():
    images = os.listdir("images")
    for im in images:
            check_size(im)
            


if __name__=="__main__":
    main()