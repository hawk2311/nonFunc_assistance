#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#define WIDTH  133
#define HEIGHT 100
#define MAX_PIXELS (WIDTH * HEIGHT)

int parse_single_int(const char *line, int *out)
{
    char *end;
    long val = strtol(line, &end, 10);

    if (line == end)
        return 0;

    while (*end) {
        if (!isspace((unsigned char)*end))
            return 0;
        end++;
    }

    *out = (int)val;
    return 1;
}

int main(void)
{
    FILE *in = fopen("output.txt", "r");
    if (!in) {
        perror("output.txt");
        return 1;
    }

    FILE *out = fopen("result.pgm", "w");
    if (!out) {
        perror("result.pgm");
        fclose(in);
        return 1;
    }

    char line[256];
    int value;
    int pixel_count = 0;

    fprintf(out, "P2\n%d %d\n255\n", WIDTH, HEIGHT);

    while (fgets(line, sizeof(line), in)) {

        /* nur reine Zahlenzeilen akzeptieren */
        if (!parse_single_int(line, &value))
            continue;

        if (pixel_count >= MAX_PIXELS)
            break;

        if (value < 0) value = 0;
        if (value > 255) value = 255;

        fprintf(out, "%d ", value);
        pixel_count++;

        if (pixel_count % WIDTH == 0)
            fprintf(out, "\n");
    }

    fclose(in);
    fclose(out);

    printf("Bild erzeugt: result.pgm (%dx%d)\n", WIDTH, HEIGHT);
    return 0;
}

