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
csv_avg = sys.argv[5] #average values of perf data


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
    with open("header/images/dataset_2/image_data_"+str(index)+".h", "w") as f:
        f.write("#ifndef IMAGE_DATA_H\n#define IMAGE_DATA_H\n\n")
        f.write("#include <stdint.h>\n\n")
        f.write(f"const uint8_t image_data[{height}][{width}] = {{\n") 
        for row in pixel_data:
            line = ", ".join(f"{val:3d}" for val in row)
            f.write(f"    {{{line}}},\n")
        f.write("};\n\n#endif\n")


def get_sizes(images):
    for im in images:
        img= Image.open("images/dataset_2/"+im).convert("L")
        image_size.append(img.size)


def create_image_data(images):
    index = 1
    for im in images:
        create_header("images/dataset_2/"+im, index)
        index+=1

def update_code(index, width, height):
    with open(name_code, "r",  encoding='utf-8') as file:
        data = file.readlines()
        data[0] = "#include \"header/images/dataset_2/image_data_"+str(index)+".h\"\n"
        data[1] = "#define WIDTH "+ str(width) +"\n"
        data[2] = "#define HEIGHT "+ str(height) + "\n"
    with open(name_code, "w",  encoding='utf-8') as file:
        file.writelines(data)
    

def collect_data(images, index):
    sum_ins = sum_cyc = sum_bra = sum_car = sum_cam = sum_time =0
    for i in range(10):
            res = subprocess.run(["perf", "stat" ,"-e" ,"instructions,cycles,branches,cache-references,cache-misses" , "./" + name_exec], capture_output=True, text=True) #executing perf stat with compiled code
            #catch relevant values
            ins = re.search("([0-9][0-9.]+)\s*(instructions)", res.stderr)
            cyc = re.search("([0-9][0-9.]+)\s*cycles", res.stderr)
            bra = re.search("([0-9][0-9.]+)\s*branches", res.stderr)
            car = re.search("([0-9][0-9.]+)\s*cache-references", res.stderr)
            cam = re.search("([0-9][0-9.]+)\s*cache-misses", res.stderr)
            time = re.search("([0-9]+,[0-9]*)\s*seconds time elapsed", res.stderr)

            #Writing to CSV File
            #with first run in whole algorithm all data written in CSV file before is overwritten
            if i<1 and index<1:
                with open(csv_data, "w") as csv_f:
                    writer = csv.writer(csv_f)
                    writer.writerow(['image_name:'+ images[index], 'image_size(width,height):'+ str(image_size[index])])
                    writer.writerow(['instructions','cycles','branches','cache-referencs','cache-misses', 'time elapsed in s'])

            #after first run all data needs to be appended
            if i<1 and index>=1:
                with open(csv_data, "a") as csv_f:
                    writer = csv.writer(csv_f)
                    writer.writerow(['image_name:'+images[index], 'image_size(width,height):'+ str(image_size[index])])
                    writer.writerow(['instructions','cycles','branches','cache-referencs','cache-misses', 'time elapsed in s'])

            with open(csv_data, "a") as csv_f:
                writer = csv.writer(csv_f)
                writer.writerow([ins.group(1), cyc.group(1), bra.group(1), car.group(1), cam.group(1), time.group(1)])

            #Creating average values
            sum_ins += int(ins.group(1).replace(".", ""))
            sum_cyc += int(cyc.group(1).replace(".", ""))
            sum_bra += int(bra.group(1).replace(".", ""))
            sum_car += int(car.group(1).replace(".", ""))
            sum_cam += int(cam.group(1).replace(".", ""))
            sum_time += float(time.group(1).replace(",", "."))


            if i >= 9:
                avg_ins = int(sum_ins/10)
                avg_cyc = int(sum_cyc/10)
                avg_bra = int(sum_bra/10)
                avg_car = int(sum_car/10)
                avg_cam = int(sum_cam/10)
                avg_time = sum_time/10

                if index<1:
                    with open(csv_avg, "w") as csv_f:
                                    writer = csv.writer(csv_f)
                                    writer.writerow(['image_name', 'image_size:(width,height)', 'instructions','cycles','branches','cache-referencs','cache-misses', 'time elapsed in s'])

                with open(csv_avg, "a") as csv_f:
                            writer = csv.writer(csv_f)
                            writer.writerow([images[index], image_size[index], avg_ins, avg_cyc, avg_bra, avg_car, avg_cam, avg_time])



    

        
    



def main():
    images = os.listdir("images/dataset_2")
    images.sort()

    #if the -header flag is set than the header files for the images will be created
    if sys.argv[6] == "-header":
         create_image_data(images)
    else:
         get_sizes(images)


    # parser = argparse.ArgumentParser()
    # parser.add_argument("-header", action="store_true")
    # args = parser.parse_args()
    # if args.header:
    #     create_image_data()
    # else:
    #     get_sizes()


    createBuffer() #with this it is possible to enter the compile command and give it to subprocess.run

    for index in range(len(image_size)):
        width, height = image_size[index] #width and height of all images got saved before
        update_code(index+1, width, height)
        subprocess.run(buffer)   #compile the code
        collect_data(images, index)
        

    


if __name__=="__main__":
    main()