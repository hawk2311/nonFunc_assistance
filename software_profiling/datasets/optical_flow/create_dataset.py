import subprocess
import os
import re
import json
import csv
import sys
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


def create_video_data(videos):
    #print(videos)
    index = 1
    for v in videos:
        video= os.path.join("./Videos", v)
        subprocess.run(["python3", "video_to_headers.py", video, str(index)])
        index+=1

def update_code(index, width, height):
    with open("optflow_2frames.cpp", "r",  encoding='utf-8') as file:
        data = file.readlines()
        data[0] = "#include \"./header/video_"+str(index)+"/image_0000.h\"\n"
        data[1] = "#include \"./header/video_"+str(index)+"/image_0001.h\"\n"
        data[2] = "#define WIDTH " + str(width)+ "\n"
        data[3] = "#define HEIGHT " + str(height) + "\n"



    with open("optflow_2frames.cpp", "w",  encoding='utf-8') as file:
        file.writelines(data)


def collect_data(videos):
    index = 1
    for v in videos:
        width, height=get_video_dimensions("./Videos/" + v)    
        update_code(index, width, height)
        subprocess.run(["g++", "-static",  "-o", "optflow_2frames" , "optflow_2frames.cpp" , "-lm"]) #compile the code 
        res = subprocess.run(["perf", "stat" ,"-e" ,"instructions,cycles,branches,cache-references,cache-misses" ,"./optflow_2frames"], capture_output=True, text=True) #executing perf stat with compiled code
        #catch relevant values
        #print(res)
        ins = re.search("([0-9][0-9.]+)\s*(instructions)", res.stderr)
        cyc = re.search("([0-9][0-9.]+)\s*cycles", res.stderr)
        bra = re.search("([0-9][0-9.]+)\s*branches", res.stderr)
        car = re.search("([0-9][0-9.]+)\s*cache-references", res.stderr)
        cam = re.search("([0-9][0-9.]+)\s*cache-misses", res.stderr)
        time = re.search("([0-9]+,[0-9]*)\s*seconds time elapsed", res.stderr)


        if index<2:
            with open("optflow_2frames_data.csv", "w") as csv_f:
                writer = csv.writer(csv_f)
                writer.writerow(['video_name', 'video_width', 'video_height', 'instructions','cycles','branches','cache-referencs','cache-misses', 'times elapsed in s'])

        with open("optflow_2frames_data.csv", "a") as csv_f:
            writer = csv.writer(csv_f)
            writer.writerow([v, width, height, ins.group(1), cyc.group(1), bra.group(1), car.group(1), cam.group(1), time.group(1)])

        index+=1

        



def main():
    #videos = os.listdir("~/home/dennis/Dokumente/Projekte/chip/Videos")
    videos = os.listdir("Videos")
    videos.sort()

    if sys.argv== "-header":
        create_video_data(videos)
        
    collect_data(videos)
        



    


if __name__=="__main__":
    main()