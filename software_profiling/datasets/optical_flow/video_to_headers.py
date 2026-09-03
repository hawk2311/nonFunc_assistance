"""
video_to_headers.py
Verwendung:
  python3 video_to_headers.py input.mp4 1

Erzeugt:
  ./header/videos/video_1/image_0000.h ...
  ./header/videos/video_1/images_index.h
"""

import os
import sys
import subprocess
import argparse
import shutil
import json
from PIL import Image
import numpy as np

BASE_OUTPUT_DIR = os.path.abspath(os.path.join("..", "header", "videos"))

def get_video_dimensions(video_in):
    cmd = [
        "ffprobe", 
        "-v", "error", 
        "-select_streams", "v:0", 
        "-show_entries", "stream=width,height", 
        "-of", "json", 
        video_in
    ]
    try:
        output = subprocess.check_output(cmd).decode("utf-8")
        data = json.loads(output)
        width = data["streams"][0]["width"]
        height = data["streams"][0]["height"]
        return width, height
    except Exception as e:
        print(f"Fehler beim Auslesen der Videogröße mit ffprobe: {e}")
        sys.exit(1)

def extract_frames(video_in, temp_frames_dir):
    os.makedirs(temp_frames_dir, exist_ok=True)
    cmd = [
        "ffmpeg", "-hide_banner", "-loglevel", "error",
        "-i", video_in,
        "-pix_fmt", "gray",
        os.path.join(temp_frames_dir, "frame%04d.png")
    ]
    #print("Extrahiere Frames aus dem Video...")
    subprocess.check_call(cmd)

def image_to_c_array(img, arr_name):
    a = np.array(img, dtype=np.uint8).flatten()
    per = 12
    lines = []
    for i in range(0, a.size, per):
        chunk = a[i:i+per]
        lines.append(", ".join(str(int(x)) for x in chunk))
    init = ",\n    ".join(lines)
    return f"const uint8_t {arr_name}[{a.size}] = {{\n    {init}\n}};\n"

def main():
    p = argparse.ArgumentParser(description="Konvertiert ein Video in C-Header-Dateien.")
    p.add_argument("video", help="Pfad zur Eingabe-Videodatei")
    p.add_argument("index", help="Index für den Zielordner (z.B. 1 für video_1)")
    args = p.parse_args()

    if not os.path.isfile(args.video):
        print(f"Fehler: Datei '{args.video}' wurde nicht gefunden.")
        sys.exit(1)

    width, height = get_video_dimensions(args.video)
    #print(f"Erkannte Videogröße: {width}x{height}")

    # Zielordner ist dynamisch basierend auf dem Index (./header/videos/video_i)
    hdr_dir = os.path.join(BASE_OUTPUT_DIR, f"video_{args.index}")
    os.makedirs(hdr_dir, exist_ok=True)

    temp_frames_dir = os.path.join(hdr_dir, "_temp_frames")

    try:
        extract_frames(args.video, temp_frames_dir)

        files = sorted([f for f in os.listdir(temp_frames_dir) if f.lower().endswith(".png")])
        if not files:
            print("Keine Frames im Video gefunden.")
            sys.exit(1)

        arr_names = []

        for idx, fname in enumerate(files):
            path = os.path.join(temp_frames_dir, fname)
            img = Image.open(path).convert("L")
            
            arr_name = f"image_{idx:04d}"
            arr_names.append(arr_name)
            
            c_code = image_to_c_array(img, arr_name)
            header_name = os.path.join(hdr_dir, f"{arr_name}.h")
            
            with open(header_name, "w") as fh:
                fh.write(f"#ifndef IMAGE_{arr_name.upper()}_H\n#define IMAGE_{arr_name.upper()}_H\n\n")
                fh.write("#include <stdint.h>\n\n")
                fh.write(c_code)
                fh.write("\n#endif\n")
            #print(f"Erstellt: {header_name}")

        idx_path = os.path.join(hdr_dir, "images_index.h")
        with open(idx_path, "w") as fh:
            fh.write("#ifndef IMAGES_INDEX_H\n#define IMAGES_INDEX_H\n\n")
            fh.write("#include <stdint.h>\n\n")
            for name in arr_names:
                fh.write(f'#include "{name}.h"\n')
            fh.write("\n")
            fh.write(f"#define IMG_WIDTH {width}\n#define IMG_HEIGHT {height}\n#define IMG_COUNT {len(arr_names)}\n\n")
            
            for name in arr_names:
                fh.write(f"extern const uint8_t {name}[{width * height}];\n")
            fh.write("\nstatic const uint8_t* const frames[IMG_COUNT] = {\n")
            for name in arr_names:
                fh.write(f"    {name},\n")
            fh.write("};\n\n#endif\n")
            
        #print(f"Erstellt: {idx_path}")
        #print(f"Fertig: {len(arr_names)} Header erstellt in {hdr_dir} (Auflösung: {width}x{height})")

    finally:
        if os.path.exists(temp_frames_dir):
            shutil.rmtree(temp_frames_dir)

if __name__ == "__main__":
    main()