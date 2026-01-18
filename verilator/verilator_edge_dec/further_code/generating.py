import numpy as np

# Eingabegröße
size = int(input("Please enter size (z.B. 600): "))
array = np.zeros((size, size), dtype=np.uint8)

# Parameter für die breite Linie
line_width = size // 5  # z.B. 1/5 der Gesamtgröße
center = size // 2      # Linie in der Mitte

half_width = line_width // 2

# Füllen des Arrays
for y in range(size):
    for x in range(size):
        dist = abs(y - center)
        if dist <= half_width:
            # Linear von 0 -> 255 zur Mitte, dann 255 -> 0
            value = int((dist / half_width) * 255)
        else:
            value = 255
        array[y, x] = value

# Optional: als C-Array exportieren
with open("sim_data.h", "w") as f:
    f.write("#ifndef SIM_DATA_H\n#define SIM_DATA_H\n\n")
    f.write("#include <stdint.h>\n\n")
    f.write(f"const uint8_t sim_data[{size}][{size}] = {{\n")
    for row in array:
        line = ", ".join(f"{val:3d}" for val in row)
        f.write(f"    {{{line}}},\n")
    f.write("};\n\n#endif\n")

print("Sim data written to pattern_data.h")
