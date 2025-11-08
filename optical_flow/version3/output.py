# #!/usr/bin/env python3
# import numpy as np
# import cv2
# import struct
# import sys
# import os

# def main():
#     if len(sys.argv) < 3:
#         print("Verwendung: python3 output.py <flow_output.raw> <output_video.mp4>")
#         sys.exit(1)

#     raw_file = sys.argv[1]
#     output_video = sys.argv[2]

#     if not os.path.exists(raw_file):
#         print(f"Fehler: Datei '{raw_file}' nicht gefunden.")
#         sys.exit(1)

#     # Datei öffnen und versuchen, Header zu lesen
#     with open(raw_file, "rb") as f:
#         header = f.read(12)
#         if len(header) < 12:
#             print("Fehler: Datei zu klein oder kein Header.")
#             sys.exit(1)

#         try:
#             w, h, n = struct.unpack("<III", header)
#         except struct.error:
#             print("Fehler: Header konnte nicht gelesen werden (falsches Format).")
#             sys.exit(1)

#         # Prüfen, ob Headerwerte plausibel sind
#         if w > 10000 or h > 10000 or n > 10000:
#             print("Warnung: Header sieht verdächtig aus. Datei könnte Text sein.")
#             f.seek(0)
#             data_text = np.loadtxt(f, dtype=np.float32)
#             total_vals = len(data_text)
#             print(f"Textdaten erkannt mit {total_vals} Werten.")
#             # Nutzer muss Dimensionen manuell anpassen
#             w = int(input("Breite w: "))
#             h = int(input("Höhe h: "))
#             n = total_vals // (2 * w * h)
#             data = data_text
#         else:
#             # Binärformat korrekt
#             data = np.frombuffer(f.read(), dtype=np.float32)

#     # Sicherstellen, dass genug Werte da sind
#     expected_vals = 2 * w * h * n
#     if len(data) < expected_vals:
#         print(f"Warnung: Erwartet {expected_vals} Floats, gefunden {len(data)} – Datei evtl. unvollständig.")
#         n = len(data) // (2 * w * h)

#     # Daten trennen in u und v
#     half = len(data) // 2
#     u = data[:half].reshape((n, h, w))
#     v = data[half:].reshape((n, h, w))

#     # Video vorbereiten
#     fourcc = cv2.VideoWriter_fourcc(*'mp4v')
#     fps = 10
#     out = cv2.VideoWriter(output_video, fourcc, fps, (w, h))

#     for i in range(n):
#         magnitude = np.sqrt(u[i] ** 2 + v[i] ** 2)
#         angle = np.arctan2(v[i], u[i])
#         hsv = np.zeros((h, w, 3), dtype=np.uint8)
#         hsv[..., 0] = ((angle + np.pi) / (2 * np.pi) * 180).astype(np.uint8)
#         hsv[..., 1] = 255
#         hsv[..., 2] = np.clip(magnitude / np.max(magnitude + 1e-6) * 255, 0, 255).astype(np.uint8)
#         rgb = cv2.cvtColor(hsv, cv2.COLOR_HSV2BGR)
#         out.write(rgb)

#     out.release()
#     print(f"Video gespeichert unter: {output_video}")

# if __name__ == "__main__":
#     main()


# # import numpy as np
# # import cv2
# # import sys
# # import os

# # def main():
# #     if len(sys.argv) < 3:
# #         print("Usage: python3 output.py flow_output.txt out.mp4")
# #         sys.exit(1)

# #     flow_file = sys.argv[1]
# #     output_file = sys.argv[2]

# #     # Datei öffnen und Zeilenweise lesen
# #     lines = []
# #     with open(flow_file, "r") as f:
# #         for line in f:
# #             # Nur Zeilen mit Zahlen akzeptieren
# #             parts = line.strip().split()
# #             if len(parts) == 2:
# #                 try:
# #                     u_val = float(parts[0])
# #                     v_val = float(parts[1])
# #                     lines.append((u_val, v_val))
# #                 except ValueError:
# #                     continue  # alles überspringen, was kein float ist

# #     if not lines:
# #         print("Keine gültigen Flow-Daten gefunden.")
# #         sys.exit(1)

# #     # Hier musst du Breite/Höhe wissen oder anpassen:
# #     WIDTH = 600
# #     HEIGHT = 600
# #     FRAME_COUNT = 1  # oder entsprechend anpassen

# #     data = np.array(lines, dtype=np.float32)
# #     if data.size < WIDTH * HEIGHT * 2:
# #         print("Warnung: Datei enthält zu wenige Werte")
# #     data = data[:WIDTH * HEIGHT, :]  # falls zu lang

# #     u = data[:, 0].reshape((HEIGHT, WIDTH))
# #     v = data[:, 1].reshape((HEIGHT, WIDTH))

# #     # Visualisierung in HSV
# #     mag, ang = cv2.cartToPolar(u, v)
# #     hsv = np.zeros((HEIGHT, WIDTH, 3), dtype=np.uint8)
# #     hsv[..., 0] = (ang * 180 / np.pi / 2).astype(np.uint8)
# #     hsv[..., 1] = 255
# #     hsv[..., 2] = np.clip(mag * 8, 0, 255).astype(np.uint8)
# #     bgr = cv2.cvtColor(hsv, cv2.COLOR_HSV2BGR)

# #     cv2.imwrite(output_file, bgr)
# #     print(f"Flow-Visualisierung in {output_file} gespeichert.")

# # if __name__ == "__main__":
# #     main()

# #!/usr/bin/env python3
# import numpy as np
# import cv2
# import re
# import sys

# if len(sys.argv) < 4:
#     print("Usage: flow_txt_to_video.py flow_output.txt width height [out.mp4]")
#     sys.exit(1)

# fname = sys.argv[1]
# W = int(sys.argv[2])
# H = int(sys.argv[3])
# outname = sys.argv[4] if len(sys.argv)>4 else "flow_vis.mp4"

# # parse text lines
# pattern = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+([-\d.eE]+)\s+([-\d.eE]+)")
# # dictionary frame -> list of tuples (x,y,u,v)
# frames = {}
# max_frame = -1
# with open(fname, "r") as f:
#     for line in f:
#         m = pattern.match(line)
#         if not m: continue
#         fr = int(m.group(1)); x=int(m.group(2)); y=int(m.group(3))
#         u=float(m.group(4)); v=float(m.group(5))
#         if fr not in frames: frames[fr] = []
#         frames[fr].append((x,y,u,v))
#         if fr > max_frame: max_frame = fr

# nframes = max_frame + 1
# print("Detected frames:", nframes)

# fourcc = cv2.VideoWriter_fourcc(*'mp4v')
# out = cv2.VideoWriter(outname, fourcc, 10, (W,H))

# for fr in range(nframes):
#     u_img = np.zeros((H,W), dtype=np.float32)
#     v_img = np.zeros((H,W), dtype=np.float32)
#     if fr in frames:
#         for (x,y,u,v) in frames[fr]:
#             if 0 <= x < W and 0 <= y < H:
#                 u_img[y,x] = u
#                 v_img[y,x] = v
#     mag, ang = cv2.cartToPolar(u_img, v_img)
#     hsv = np.zeros((H,W,3), dtype=np.uint8)
#     hsv[...,0] = (ang * 180.0 / np.pi / 2.0).astype(np.uint8)
#     hsv[...,1] = 255
#     mmax = np.max(mag) if np.max(mag) > 1e-6 else 1.0
#     hsv[...,2] = (np.clip(mag / mmax * 255.0, 0, 255)).astype(np.uint8)
#     bgr = cv2.cvtColor(hsv, cv2.COLOR_HSV2BGR)
#     out.write(bgr)

# out.release()
# print("Wrote", outname)


import numpy as np
import cv2

# === KONFIGURATION ===
WIDTH  = 600
HEIGHT = 600
OUTPUT_VIDEO = "out.mp4"
INPUT_FILE   = "flow_output.txt"
FPS = 10  # Frames pro Sekunde

def flow_to_color(u, v):
    """
    Wandelt optischen Flow (u,v) in ein Farbbild um.
    Hue = Richtung, Value = Stärke.
    """
    h, w = u.shape
    mag, ang = cv2.cartToPolar(u, v)
    hsv = np.zeros((h, w, 3), dtype=np.uint8)
    hsv[..., 0] = ((ang + np.pi) / (2 * np.pi) * 180).astype(np.uint8)  # Richtung → Farbe
    hsv[..., 1] = 255  # volle Sättigung
    hsv[..., 2] = cv2.normalize(mag, None, 0, 255, cv2.NORM_MINMAX)  # Stärke → Helligkeit
    return cv2.cvtColor(hsv, cv2.COLOR_HSV2BGR)


def parse_flow_file(filename):
    """Liest flow_output.txt und erzeugt Liste von (u,v)-Arrays."""
    flows = []
    with open(filename, "r") as f:
        lines = f.readlines()

    current_u = []
    current_v = []
    reading = False

    for line in lines:
        if line.startswith("START_FLOW"):
            reading = True
            current_u = []
            current_v = []
        elif line.startswith("END_FLOW"):
            if current_u and current_v:
                size = len(current_u)
                dim = int(np.sqrt(size))
                if dim * dim != size:
                    print(f"  Warnung: Ungerades Flow-Feld ({size} Werte, nicht quadratisch).")
                u_arr = np.array(current_u, dtype=np.float32).reshape((dim, dim))
                v_arr = np.array(current_v, dtype=np.float32).reshape((dim, dim))
                flows.append((u_arr, v_arr))
            reading = False
        elif reading:
            try:
                u_val, v_val = map(float, line.strip().split())
                current_u.append(u_val)
                current_v.append(v_val)
            except ValueError:
                continue  # Zeilen überspringen, die keine Zahlen enthalten

    return flows


def main():
    flows = parse_flow_file(INPUT_FILE)
    print(f"{len(flows)} Flow-Felder gefunden.")

    if not flows:
        print("Keine Flow-Daten gefunden – prüfe die Eingabedatei!")
        return

    # Video Writer
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    out = cv2.VideoWriter(OUTPUT_VIDEO, fourcc, FPS, (WIDTH, HEIGHT))

    for i, (u, v) in enumerate(flows):
        frame = flow_to_color(u, v)
        out.write(frame)
        print(f"Frame {i+1}/{len(flows)} verarbeitet")

    out.release()
    print(f"Video gespeichert als: {OUTPUT_VIDEO}")

if __name__ == "__main__":
    main()
