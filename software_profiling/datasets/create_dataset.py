import subprocess
import os
from PIL import Image
import numpy as np

image_size = [] #store size of images for later use


def create_header(name, index):
    img = Image.open(name).convert("L")
    width, height = img.size # returns #tupel(width, height)
    image_size.append(img.size)
    #convert into NumPy-Array 
    pixel_data = np.array(img, dtype=np.uint8) 
    # export it as a C-array
    with open("header/image_data_"+str(index)+".h", "w") as f:
        f.write("#ifndef IMAGE_DATA_H\n#define IMAGE_DATA_H\n\n")
        f.write("#include <stdint.h>\n\n")
        f.write(f"const uint8_t image_data[{height}][{width}] = {{\n") #600 was 64 before changing
        for row in pixel_data:
            line = ", ".join(f"{val:3d}" for val in row)
            f.write(f"    {{{line}}},\n")
        f.write("};\n\n#endif\n")

def create_image_data():
    images = os.listdir("images")
    images.sort()
    index = 1
    for im in images:
        create_header("images/"+im, index)
        index+=1

def update_code(index, width, height):
    with open("edgeBare.c", "r",  encoding='utf-8') as file:
        data = file.readlines()
        data[0] = "#include image_harder_"+str(index)+".h\n"
        data[1] = "#define WIDTH "+ str(width) +"\n"
        data[2] = "#define HEIGHT "+ str(height) + "\n"
    with open("edgeBare.c", "w",  encoding='utf-8') as file:
        file.writelines(data)

def collect_data():
    index = 0
    #for index in len(image_size):
    width, height = image_size[index]
    update_code(index+1, width, height)
    #TODO: compile code, run code with perf, write in CSV, do next one
    



#Ausführung von Befehl
#subprocess.run(["ls"])

def main():
    create_image_data()
    collect_data()
        



    


if __name__=="__main__":
    main()