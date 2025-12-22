#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "sim_data.h"

const int width = 100;
const int height = 100;



int main(){


    uint8_t gray[width*height];

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            gray[y * width + x] = sim_data[x][y];
        }
    }

    if (stbi_write_png("edges.png", width, height, 1, gray, width)) {
        printf("result saved in: edges.png\n");
    } else {
        printf("Warning: Result could not be saved .\n");
    }


return 0;

}