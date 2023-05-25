#include <stdint.h>
#include <stdbool.h>

// #define SCREEN_WIDTH 640
// #define SCREEN_HEIGHT 480

// #define SCREEN_WIDTH 256
// #define SCREEN_HEIGHT 256

#define SCREEN_FACTOR 2
#define SCREEN_WIDTH  (64*SCREEN_FACTOR)
#define SCREEN_HEIGHT (64*SCREEN_FACTOR)

// #define SCREEN_WIDTH 2560
// #define SCREEN_HEIGHT 1600

void init_sdl();
void update_sdl(uint32_t *pixels);
void close_sdl();

void init_x264(char * filename);
int update_x264(uint32_t *pixels);
void close_x264();
