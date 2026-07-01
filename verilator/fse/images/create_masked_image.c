#include <stdio.h>
#include <stdlib.h>
#include "image_data.h"
#include "mask_data.h"
//#include "mask_image.h"

#define SIZE 200

void write_header_2d(int data[SIZE][SIZE], const char *filename) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        perror("fopen");
        return;
    }

    fprintf(f, "#ifndef MASK_IMAGE_H\n#define MASK_IMAGE_H\n\n");
    fprintf(f, "static int mask_image[%d][%d] = {\n", SIZE, SIZE);

    for (int i = 0; i < SIZE; i++) {
        fprintf(f, "    {");
        for (int j = 0; j < SIZE; j++) {
            fprintf(f, "%d%s", data[i][j], (j < SIZE - 1) ? "," : "");
        }
        fprintf(f, "}%s\n", (i < SIZE - 1) ? "," : "");
    }

    fprintf(f, "};\n\n#endif // MASK_IMAGE_H\n");
    fclose(f);
}

int main(void) {

    int mask_image[SIZE][SIZE];

    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            if (mask_data[i][j] <= 10) {
                mask_image[i][j] = mask_data[i][j];
            } else {
                mask_image[i][j] = image_data[i][j];
            }
        }
    }

    write_header_2d(mask_image, "mask_image.h");

    return 0;
}
