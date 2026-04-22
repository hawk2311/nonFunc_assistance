#!/usr/bin/env python3
import cv2
import numpy as np
import sys
import os

def main():
    # --- Argumente prüfen ---
    if len(sys.argv) < 3:
        print("Verwendung: python3 video_to_h.py <input_video> <output_h>")
        sys.exit(1)

    video_in = sys.argv[1]
    h_out    = sys.argv[2]

    if not os.path.exists(video_in):
        print(f"Fehler: Eingabedatei '{video_in}' nicht gefunden.")
        sys.exit(1)

    # --- Video öffnen ---
    cap = cv2.VideoCapture(video_in)
    if not cap.isOpened():
        print(f"Fehler: Konnte '{video_in}' nicht öffnen.")
        sys.exit(1)

    # --- Metadaten ---
    width  = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    frame_count_meta = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))

    print(f"Metadaten: {width}x{height}, Frames (Meta): {frame_count_meta}")

    # --- Frames lesen ---
    frames = []
    while True:
        ret, frame = cap.read()
        if not ret:
            break

        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        frames.append(gray)

    cap.release()

    if len(frames) == 0:
        print("Fehler: Keine Frames im Video gefunden.")
        sys.exit(1)

    frames = np.array(frames, dtype=np.uint8)
    num_frames, h, w = frames.shape

    print(f"Video gelesen: {num_frames} Frames, Auflösung: {w}x{h}")

    # --- .h-Datei schreiben ---
    array_name = os.path.splitext(os.path.basename(video_in))[0]
    array_name = array_name.replace("-", "_")

    with open(h_out, "w") as f:
        f.write("#pragma once\n\n")
        f.write("#include <stdint.h>\n\n")

        f.write(f"#define VIDEO_T {num_frames}\n")
        f.write(f"#define VIDEO_H {h}\n")
        f.write(f"#define VIDEO_W {w}\n\n")

        f.write(f"const uint8_t video_data"
                f"[VIDEO_T][VIDEO_H][VIDEO_W] = {{\n")

        for t in range(num_frames):
            f.write("  {\n")
            for y in range(h):
                row = ", ".join(str(int(v)) for v in frames[t][y])
                f.write(f"    {{{row}}},\n")
            f.write("  },\n")

        f.write("};\n")

    print(f"C-Header erzeugt: {h_out}")

if __name__ == "__main__":
    main()
