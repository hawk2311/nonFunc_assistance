
import numpy as np
from PIL import Image

#WIDTH = 1280
#HEIGHT = 720
WIDTH = 426
HEIGHT = 240

data = np.loadtxt("flow_output.txt", dtype=np.uint8)

img = data.reshape((HEIGHT, WIDTH, 3))

print("Enter name for image file:")

name = input()
Image.fromarray(img).save(name)
