import glob
import cv2

files = sorted(glob.glob("flow_vis_*.ppm"))

frames = []
for f in files:
    img = cv2.imread(f)
    frames.append(img)

h, w, _ = frames[0].shape

out = cv2.VideoWriter(
    "flow_video.mp4",
    cv2.VideoWriter_fourcc(*'mp4v'),
    10,
    (w,h)
)

for f in frames:
    out.write(f)

out.release()
