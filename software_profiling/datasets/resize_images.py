import subprocess
import os
import re
import csv
import sys
from PIL import Image
import numpy as np

image_size=[]

def get_sizes(images):
    for im in images:
        img= Image.open("images/"+im).convert("L")
        image_size.append(img.size)


def main():
    images = os.listdir("images")
    for im in images:
            img= Image.open("images/"+im).convert("L")
            #image_size.append(img.size)
            if img.size[0] > 700 or img.size[0] > 700:
                print(im,":",img.size)

    print(image_size)

if __name__=="__main__":
    main()