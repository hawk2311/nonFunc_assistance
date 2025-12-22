# flow_to_video.py
import numpy as np, cv2

flow = np.load("flow_output.npy")  # (frames, h, w, 2)
h, w = flow.shape[1:3]
out = cv2.VideoWriter("flow_vis.mp4", cv2.VideoWriter_fourcc(*'mp4v'), 25, (w, h))

for f in flow:
    mag, ang = cv2.cartToPolar(f[...,0], f[...,1])
    hsv = np.zeros((h,w,3), np.uint8)
    hsv[...,0] = ang * 180/np.pi / 2
    hsv[...,1] = 255
    hsv[...,2] = cv2.normalize(mag, None, 0, 255, cv2.NORM_MINMAX)
    rgb = cv2.cvtColor(hsv, cv2.COLOR_HSV2BGR)
    out.write(rgb)

out.release()
