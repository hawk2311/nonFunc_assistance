import cv2
import numpy as np

# Bilder laden
img1 = cv2.imread("input_videos/split/split_video/frames/frame0001.png", cv2.IMREAD_GRAYSCALE)
img2 = cv2.imread("input_videos/split/split_video/frames/frame0002.png", cv2.IMREAD_GRAYSCALE)

if img1 is None or img2 is None:
    print("Error loading images")
    exit()

# Features finden
p0 = cv2.goodFeaturesToTrack(img1, maxCorners=500, qualityLevel=0.01, minDistance=10)

if p0 is None:
    print("No features found")
    exit()

# Lucas-Kanade Optical Flow
p1, status, err = cv2.calcOpticalFlowPyrLK(img1, img2, p0, None)

# Nur gültige Punkte
good_old = p0[status.flatten() == 1]
good_new = p1[status.flatten() == 1]

# reshape in (N,2)
good_old = good_old.reshape(-1, 2)
good_new = good_new.reshape(-1, 2)

# Bewegung berechnen
flow = good_new - good_old
vx = flow[:, 0]
vy = flow[:, 1]
mag = np.sqrt(vx**2 + vy**2)

# Durchschnitt über alle validen Punkte
avg_vx = np.mean(vx)
avg_vy = np.mean(vy)
avg_mag = np.mean(mag)

print("\n[RESULT] Lucas-Kanade flow:")
print("vx =", avg_vx)
print("vy =", avg_vy)
print("magnitude =", avg_mag)