#pragma once

#include <stdint.h>

#define VIDEO_T 17
#define VIDEO_H 720
#define VIDEO_W 1280

uint8_t optical_flow_data[VIDEO_T-1][VIDEO_H][VIDEO_W] = {};
