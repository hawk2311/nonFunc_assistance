#!/usr/bin/env python3
"""
video_to_frames_and_headers.py
Usage:
  python3 video_to_frames_and_headers.py input.mp4 out_dir --width 128 --height 128

Output:
  out_dir/frame0000.png    (optional)
  out_dir/image_0000.h
  ...
  out_dir/images_index.h   (includes all image headers and defines FRAME_COUNT, WIDTH, HEIGHT)
"""
import os, sys, subprocess, argparse
from PIL import Image
import numpy as np

def extract_frames(video_in, frames_dir):
    os.makedirs(frames_dir, exist_ok=True)
    # Use ffmpeg to extract grayscale frames as PNGs
    cmd = [
        "ffmpeg", "-hide_banner", "-loglevel", "error",
        "-i", video_in,
        "-vf", "fps=25", # optional: change framerate if needed
        "-pix_fmt", "gray",
        os.path.join(frames_dir, "frame%04d.png")
    ]
    print("Running:", " ".join(cmd))
    subprocess.check_call(cmd)

def image_to_c_array(img, arr_name):
    # img: PIL grayscale image (mode 'L'), already resized to W x H
    a = np.array(img, dtype=np.uint8).flatten()
    # produce C initializer
    lines = []
    per_line = 16
    for i in range(0, a.size, per_line):
        chunk = a[i:i+per_line]
        lines.append(", ".join(str(int(x)) for x in chunk))
    init = ",\n    ".join(lines)
    return f"const uint8_t {arr_name}[{a.size}] = {{\n    {init}\n}};\n"

def main():
    p = argparse.ArgumentParser()
    p.add_argument("video", help="input video")
    p.add_argument("outdir", help="output dir for frames and headers")
    p.add_argument("--width", type=int, default=128)
    p.add_argument("--height", type=int, default=128)
    p.add_argument("--no-extract", action="store_true", help="skip ffmpeg extract (already have png frames)")
    args = p.parse_args()

    frames_dir = os.path.join(args.outdir, "frames")
    os.makedirs(args.outdir, exist_ok=True)
    if not args.no_extract:
        extract_frames(args.video, frames_dir)
    else:
        if not os.path.isdir(frames_dir):
            raise SystemExit("frames dir missing; use --no-extract only if frames exist")

    # gather PNG files sorted
    files = sorted([f for f in os.listdir(frames_dir) if f.lower().endswith(".png")])
    if len(files) < 2:
        raise SystemExit("Need at least two frames")

    hdr_dir = os.path.join(args.outdir, "headers")
    os.makedirs(hdr_dir, exist_ok=True)

    arr_names = []
    for idx, fname in enumerate(files):
        img_path = os.path.join(frames_dir, fname)
        img = Image.open(img_path).convert("L")
        img = img.resize((args.width, args.height), Image.BILINEAR)
        arr_name = f"image_{idx:04d}"
        arr_names.append(arr_name)
        c_code = image_to_c_array(img, arr_name)
        header_name = os.path.join(hdr_dir, f"{arr_name}.h")
        with open(header_name, "w") as fh:
            fh.write("#ifndef IMAGE_%s_H\n#define IMAGE_%s_H\n\n" % (arr_name.upper(), arr_name.upper()))
            fh.write("#include <stdint.h>\n\n")
            fh.write(c_code)
            fh.write("\n#endif\n")
        print("Wrote", header_name)

    # make images_index.h which includes all headers and exposes frames array
    idx_path = os.path.join(hdr_dir, "images_index.h")
    with open(idx_path, "w") as fh:
        fh.write("#ifndef IMAGES_INDEX_H\n#define IMAGES_INDEX_H\n\n")
        fh.write("#include <stdint.h>\n\n")
        for name in arr_names:
            fh.write('#include "%s.h"\n' % name)
        fh.write("\n")
        W = args.width
        H = args.height
        N = len(arr_names)
        fh.write(f"#define IMG_WIDTH {W}\n#define IMG_HEIGHT {H}\n#define IMG_COUNT {N}\n\n")
        # declare extern arrays (they are const defined in their headers)
        for name in arr_names:
            fh.write(f"extern const uint8_t {name}[{W*H}];\n")
        fh.write("\n/* frames as pointers to arrays (each is a uint8_t pointer to W*H values) */\n")
        fh.write(f"static const uint8_t* const frames[{N}] = {{\n")
        for name in arr_names:
            fh.write(f"    {name},\n")
        fh.write("};\n\n")
        fh.write("#endif\n")
    print("Wrote", idx_path)
    print("Done. Created %d headers, width=%d height=%d" % (N, W, H))

if __name__ == "__main__":
    main()
