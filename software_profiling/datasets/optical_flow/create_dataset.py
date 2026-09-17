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


def create_video_data(videos):
    #print(videos)
    index = 1
    for v in videos:
        video= os.path.join("../videos", v)
        subprocess.run(["python3", "video_to_headers.py", video, str(index)])
        index+=1

def update_code(index):
    with open("optflow.cpp", "r",  encoding='utf-8') as file:
        data = file.readlines()
        data[0] = "#include \"../header/videos/video_"+str(index)+"/images_index.h\"\n"
    with open("optflow.cpp", "w",  encoding='utf-8') as file:
        file.writelines(data)

#TODO
def collect_data(videos):
    index = 1
    for v in videos:    
        update_code(index)
        subprocess.run(["g++", "-static",  "-o", "optflow" , "optflow.cpp" , "-lm"]) #compile the code 
        res = subprocess.run(["perf", "stat" ,"-e" ,"instructions,cycles,branches,cache-references,cache-misses" ,"./optflow"], capture_output=True, text=True) #executing perf stat with compiled code
        #catch relevant values
        #print(res)
        ins = re.search("([0-9][0-9.]+)\s*(instructions)", res.stderr)
        cyc = re.search("([0-9][0-9.]+)\s*cycles", res.stderr)
        bra = re.search("([0-9][0-9.]+)\s*branches", res.stderr)
        car = re.search("([0-9][0-9.]+)\s*cache-references", res.stderr)
        cam = re.search("([0-9][0-9.]+)\s*cache-misses", res.stderr)

        if index<2:
            with open("optflow_data.csv", "w") as csv_f:
                writer = csv.writer(csv_f)
                writer.writerow(['video_name','instructions','cycles','branches','cache-referencs','cache-misses'])

        with open("optflow_data.csv", "a") as csv_f:
            writer = csv.writer(csv_f)
            writer.writerow([v, ins.group(1), cyc.group(1), bra.group(1), car.group(1), cam.group(1)])

        index+=1

        



def main():
    videos = os.listdir("../videos")
    videos.sort()


    #create_video_data(videos)
    collect_data(videos)
        



    


if __name__=="__main__":
    main()