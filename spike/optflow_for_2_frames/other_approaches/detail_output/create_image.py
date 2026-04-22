import numpy as np
from PIL import Image

# Parameter (müssen exakt zu deinem C-Code passen!)
#WIDTH = 1280
#HEIGHT = 720
WIDTH = 426
HEIGHT = 240

# 1. Datei einlesen
with open("flow_output.txt", "r") as f:
    data = f.readlines()

# 2. In Integer umwandeln
pixels = np.array([int(x.strip()) for x in data], dtype=np.uint8)

# 3. In Bildform bringen
image = pixels.reshape((HEIGHT, WIDTH))

# 4. Bild speichern
img = Image.fromarray(image, mode='L')  # 'L' = Graustufen
img.save("flow_from_txt.png")

print("Bild gespeichert als flow_from_txt_2.png")
