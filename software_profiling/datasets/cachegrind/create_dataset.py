import subprocess
import os
import re
import csv
import sys
from PIL import Image
import numpy as np


image_size = [] #store size of images for later use

comp_cmd=sys.argv[1] #you need to add the compilation command for your code 
name_code = sys.argv[2]
name_exec = sys.argv[3] #add the name of the executable of the compilation, must be the same as in the command before
csv_data = sys.argv[4] #output for perf data


buffer= []
a = 0
def addNextWort():
    global a
    temp = []
    # stop with end of command
    while a < len(comp_cmd) and comp_cmd[a] != ' ' :
        temp.append(comp_cmd[a])
        a += 1
    # ignore spaces or dots
    while a < len(comp_cmd) and comp_cmd[a] == ' ':
        a += 1
    return ''.join(temp)

def createBuffer():
    global a
    while a < len(comp_cmd):
        next_wort = addNextWort()
        if next_wort:  # ignore empty strings
            buffer.append(next_wort)

    return buffer

def create_header(name, index):
    img = Image.open(name).convert("L")
    width, height = img.size # returns #tupel(width, height)
    image_size.append(img.size)
    #convert into NumPy-Array 
    pixel_data = np.array(img, dtype=np.uint8) 
    # export it as a C-array
    with open("../header/images/dataset_2/image_data_"+str(index)+".h", "w") as f:
        f.write("#ifndef IMAGE_DATA_H\n#define IMAGE_DATA_H\n\n")
        f.write("#include <stdint.h>\n\n")
        f.write(f"const uint8_t image_data[{height}][{width}] = {{\n") 
        for row in pixel_data:
            line = ", ".join(f"{val:3d}" for val in row)
            f.write(f"    {{{line}}},\n")
        f.write("};\n\n#endif\n")


def get_sizes(images):
    for im in images:
        img= Image.open("../images/dataset_2/"+im).convert("L")
        image_size.append(img.size)


def create_image_data(images):
    index = 1
    for im in images:
        create_header("../images/dataset_2/"+im, index)
        index+=1

def update_code(index, width, height):
    with open(name_code, "r",  encoding='utf-8') as file:
        data = file.readlines()
        data[0] = "#include \"../header/images/dataset_2/image_data_"+str(index)+".h\"\n"
        data[1] = "#define WIDTH "+ str(width) +"\n"
        data[2] = "#define HEIGHT "+ str(height) + "\n"
    with open(name_code, "w",  encoding='utf-8') as file:
        file.writelines(data)
    

def collect_data(images, index):
    res = subprocess.run(["valgrind", "--tool=cachegrind" ,"--branch-sim=yes" , "./" + name_exec], capture_output=True, text=True) #executing perf stat with compiled code
    #catch relevant values
    #print(res.stderr)
    ins = re.search("(I\srefs:)\s*([0-9,]+)", res.stderr)
    bra = re.search("(Branches:)\s*([0-9,]+)", res.stderr)

    #Writing to CSV File
    #with first run in whole algorithm all data written in CSV file before is overwritten
    if index<1:
        with open(csv_data, "w") as csv_f:
            writer = csv.writer(csv_f)
            writer.writerow(['image_name:'+ images[index], 'image_size(width,height):'+ str(image_size[index])])
            writer.writerow(['instructions','branches'])

    #after first run all data needs to be appended
    if index>=1:
        with open(csv_data, "a") as csv_f:
            writer = csv.writer(csv_f)
            writer.writerow(['image_name:'+images[index], 'image_size(width,height):'+ str(image_size[index])])
            writer.writerow(['instructions','branches'])

    with open(csv_data, "a") as csv_f:
        writer = csv.writer(csv_f)
        writer.writerow([ins.group(2), bra.group(2)])

         


    

        
    



def main():
    images = os.listdir("../images/dataset_2")
    images.sort()

    #TODO
    #if the -header flag is set than the header files for the images will be created
    # if sys.argv[6] == "-header":
    #      create_image_data(images)
    # else:
    #      get_sizes(images)


    # parser = argparse.ArgumentParser()
    # parser.add_argument("-header", action="store_true")
    # args = parser.parse_args()
    # if args.header:
    #     create_image_data()
    # else:
    #     get_sizes()

    create_image_data(images)
    createBuffer() #with this it is possible to enter the compile command and give it to subprocess.run

    for index in range(len(image_size)):
        width, height = image_size[index] #width and height of all images got saved before
        update_code(index+1, width, height)
        subprocess.run(buffer)   #compile the code
        collect_data(images, index)
        

    subprocess.run(["rm", "cachegrind.out.*"])



if __name__=="__main__":
    main()