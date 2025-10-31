#!/usr/bin/env python3
import numpy as np
import cv2
import struct
import sys
import os

def main():
    if len(sys.argv) < 3:
        print("Verwendung: python3 output.py <flow_output.raw> <output_video.mp4>")
        sys.exit(1)

    raw_file = sys.argv[1]
    output_video = sys.argv[2]

    if not os.path.exists(raw_file):
        print(f"Fehler: Datei '{raw_file}' nicht gefunden.")
        sys.exit(1)

    # Datei öffnen und versuchen, Header zu lesen
    with open(raw_file, "rb") as f:
        header = f.read(12)
        if len(header) < 12:
            print("Fehler: Datei zu klein oder kein Header.")
            sys.exit(1)

        try:
            w, h, n = struct.unpack("<III", header)
        except struct.error:
            print("Fehler: Header konnte nicht gelesen werden (falsches Format).")
            sys.exit(1)

        # Prüfen, ob Headerwerte plausibel sind
        if w > 10000 or h > 10000 or n > 10000:
            print("Warnung: Header sieht verdächtig aus. Datei könnte Text sein.")
            f.seek(0)
            data_text = np.loadtxt(f, dtype=np.float32)
            total_vals = len(data_text)
            print(f"Textdaten erkannt mit {total_vals} Werten.")
            # Nutzer muss Dimensionen manuell anpassen
            w = int(input("Breite w: "))
            h = int(input("Höhe h: "))
            n = total_vals // (2 * w * h)
            data = data_text
        else:
            # Binärformat korrekt
            data = np.frombuffer(f.read(), dtype=np.float32)

    # Sicherstellen, dass genug Werte da sind
    expected_vals = 2 * w * h * n
    if len(data) < expected_vals:
        print(f"Warnung: Erwartet {expected_vals} Floats, gefunden {len(data)} – Datei evtl. unvollständig.")
        n = len(data) // (2 * w * h)

    # Daten trennen in u und v
    half = len(data) // 2
    u = data[:half].reshape((n, h, w))
    v = data[half:].reshape((n, h, w))

    # Video vorbereiten
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    fps = 10
    out = cv2.VideoWriter(output_video, fourcc, fps, (w, h))

    for i in range(n):
        magnitude = np.sqrt(u[i] ** 2 + v[i] ** 2)
        angle = np.arctan2(v[i], u[i])
        hsv = np.zeros((h, w, 3), dtype=np.uint8)
        hsv[..., 0] = ((angle + np.pi) / (2 * np.pi) * 180).astype(np.uint8)
        hsv[..., 1] = 255
        hsv[..., 2] = np.clip(magnitude / np.max(magnitude + 1e-6) * 255, 0, 255).astype(np.uint8)
        rgb = cv2.cvtColor(hsv, cv2.COLOR_HSV2BGR)
        out.write(rgb)

    out.release()
    print(f"Video gespeichert unter: {output_video}")

if __name__ == "__main__":
    main()


# import numpy as np
# import cv2
# import sys
# import os

# def main():
#     if len(sys.argv) < 3:
#         print("Usage: python3 output.py flow_output.txt out.mp4")
#         sys.exit(1)

#     flow_file = sys.argv[1]
#     output_file = sys.argv[2]

#     # Datei öffnen und Zeilenweise lesen
#     lines = []
#     with open(flow_file, "r") as f:
#         for line in f:
#             # Nur Zeilen mit Zahlen akzeptieren
#             parts = line.strip().split()
#             if len(parts) == 2:
#                 try:
#                     u_val = float(parts[0])
#                     v_val = float(parts[1])
#                     lines.append((u_val, v_val))
#                 except ValueError:
#                     continue  # alles überspringen, was kein float ist

#     if not lines:
#         print("Keine gültigen Flow-Daten gefunden.")
#         sys.exit(1)

#     # Hier musst du Breite/Höhe wissen oder anpassen:
#     WIDTH = 600
#     HEIGHT = 600
#     FRAME_COUNT = 1  # oder entsprechend anpassen

#     data = np.array(lines, dtype=np.float32)
#     if data.size < WIDTH * HEIGHT * 2:
#         print("Warnung: Datei enthält zu wenige Werte")
#     data = data[:WIDTH * HEIGHT, :]  # falls zu lang

#     u = data[:, 0].reshape((HEIGHT, WIDTH))
#     v = data[:, 1].reshape((HEIGHT, WIDTH))

#     # Visualisierung in HSV
#     mag, ang = cv2.cartToPolar(u, v)
#     hsv = np.zeros((HEIGHT, WIDTH, 3), dtype=np.uint8)
#     hsv[..., 0] = (ang * 180 / np.pi / 2).astype(np.uint8)
#     hsv[..., 1] = 255
#     hsv[..., 2] = np.clip(mag * 8, 0, 255).astype(np.uint8)
#     bgr = cv2.cvtColor(hsv, cv2.COLOR_HSV2BGR)

#     cv2.imwrite(output_file, bgr)
#     print(f"Flow-Visualisierung in {output_file} gespeichert.")

# if __name__ == "__main__":
#     main()
