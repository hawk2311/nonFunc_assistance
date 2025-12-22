#!/usr/bin/env python3
import cv2
import numpy as np
import struct
import sys
import os

def main():
    # --- Argumente prüfen ---
    if len(sys.argv) < 3:
        print("Verwendung: python3 video_to_raw.py <input_video> <output_raw>")
        sys.exit(1)

    video_in = sys.argv[1]
    raw_out = sys.argv[2]

    if not os.path.exists(video_in):
        print(f"Fehler: Eingabedatei '{video_in}' nicht gefunden.")
        sys.exit(1)

    # --- Video öffnen ---
    cap = cv2.VideoCapture(video_in)
    if not cap.isOpened():
        print(f"Fehler: Konnte '{video_in}' nicht öffnen.")
        sys.exit(1)

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

    # --- Header + Frames in .raw-Datei schreiben ---
    with open(raw_out, "wb") as f:
        # Header: width, height, num_frames (je 4 Bytes, Little Endian)
        f.write(struct.pack("<III", w, h, num_frames))
        # Daten: Grauwertframes
        f.write(frames.tobytes())

    print(f"Rohdaten gespeichert in: {raw_out}")

if __name__ == "__main__":
    main()
