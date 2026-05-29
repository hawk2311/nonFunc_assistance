This directory shows the idea of splitting a video into single frames. This is currently used in this repo in Optical Flow. It is possible to generate the frames of a video with the help of the image.py file.
Use this file with the command:
python3 video_to_headers.py input.mp4 outdir --width 128 --height 128 --max-frames 50

This will generate outdir/frames where you can find all the frames of the given video and outdir/headers which contain the data of the single frames.

Edited: 29.05.2026
