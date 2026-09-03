import subprocess
import os
import re
import json
import csv
from PIL import Image
import numpy as np

def get_video_dimensions(video_in):
    cmd = [
        "ffprobe", 
        "-v", "error", 
        "-select_streams", "v:0", 
        "-show_entries", "stream=width,height", 
        "-of", "json", 
        video_in
    ]
    
    output = subprocess.check_output(cmd).decode("utf-8")
    data = json.loads(output)
    width = data["streams"][0]["width"]
    height = data["streams"][0]["height"]
    return width, height


def create_video_data():
    videos = os.listdir("../videos")
    videos.sort()
    #print(videos)
    index = 1
    for v in videos:
        video= os.path.join("../videos", v)
        subprocess.run(["python3", "video_to_headers.py", video, str(index)])
        index+=1

#TODO
def update_code(index, width, height):
    with open("fftw.c", "r",  encoding='utf-8') as file:
        data = file.readlines()
        data[0] = "#include \"header/image_data_"+str(index)+".h\"\n"
        data[1] = "#define WIDTH "+ str(width) +"\n"
        data[2] = "#define HEIGHT "+ str(height) + "\n"
    with open("fftw.c", "w",  encoding='utf-8') as file:
        file.writelines(data)

#TODO
def collect_data():
    #index = 0
   
        subprocess.run(["gcc", "-static",  "-o", "fftw" , "fftw.c" , "-lfftw3", "-lm"]) #compile the code
        #riscv64-unknown-elf-gcc -static   -I./fftw-3.3.10/install/include   -L./fftw-3.3.10/install/lib -o fftw.elf fftw.c -lfftw3 -lm
        res = subprocess.run(["perf", "stat" ,"-e" ,"instructions,cycles,branches,cache-references,cache-misses" ,"./fftw"], capture_output=True, text=True) #executing perf stat with compiled code
        #catch relevant values
        ins = re.search("([0-9][0-9.]+)\s*(instructions)", res.stderr)
        cyc = re.search("([0-9][0-9.]+)\s*cycles", res.stderr)
        bra = re.search("([0-9][0-9.]+)\s*branches", res.stderr)
        car = re.search("([0-9][0-9.]+)\s*cache-references", res.stderr)
        cam = re.search("([0-9][0-9.]+)\s*cache-misses", res.stderr)

        # if index<1:
        #     with open("fftw_data.csv", "w") as csv_f:
        #         writer = csv.writer(csv_f)
        #         writer.writerow(['image_name','instructions','cycles','branches','cache-referencs','cache-misses'])

        # with open("fftw_data.csv", "a") as csv_f:
        #     writer = csv.writer(csv_f)
        #     writer.writerow(["kodim"+str(index+1), ins.group(1), cyc.group(1), bra.group(1), car.group(1), cam.group(1)])





def main():
    create_video_data()
    #collect_data()
        



    


if __name__=="__main__":
    main()