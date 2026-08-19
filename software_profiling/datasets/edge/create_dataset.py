import subprocess
import os
import re
import csv
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
        f.write(f"const uint8_t image_data[{height}][{width}] = {{\n") 
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
        data[0] = "#include \"header/image_data_"+str(index)+".h\"\n"
        data[1] = "#define WIDTH "+ str(width) +"\n"
        data[2] = "#define HEIGHT "+ str(height) + "\n"
    with open("edgeBare.c", "w",  encoding='utf-8') as file:
        file.writelines(data)

def collect_data():
    #index = 0
    for index in range(len(image_size)):
        width, height = image_size[index] #width and height of all images got saved before
        update_code(index+1, width, height) 
        subprocess.run(["gcc", "-static",  "-o", "edgeBare" , "edgeBare.c" , "-lm"]) #compile the code
        res = subprocess.run(["perf", "stat" ,"-e" ,"instructions,cycles,branches,cache-references,cache-misses" ,"./edgeBare"], capture_output=True, text=True) #executing perf stat with compiled code
        #catch relevant values
        ins = re.search("([0-9][0-9.]+)\s*(instructions)", res.stderr)
        cyc = re.search("([0-9][0-9.]+)\s*cycles", res.stderr)
        bra = re.search("([0-9][0-9.]+)\s*branches", res.stderr)
        car = re.search("([0-9][0-9.]+)\s*cache-references", res.stderr)
        cam = re.search("([0-9][0-9.]+)\s*cache-misses", res.stderr)

        if index<1:
            with open("edge_dec_data.csv", "w") as csv_f:
                writer = csv.writer(csv_f)
                writer.writerow(['image_name','instructions','cycles','branches','cache-referencs','cache-misses'])

        with open("edge_dec_data.csv", "a") as csv_f:
            writer = csv.writer(csv_f)
            writer.writerow(["kodim"+str(index+1), ins.group(1), cyc.group(1), bra.group(1), car.group(1), cam.group(1)])





def main():
    create_image_data()
    collect_data()
        



    


if __name__=="__main__":
    main()