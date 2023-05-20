#include <stdio.h>
#include <stdlib.h>

#include "sim/sim.h"
#include "const.h"

// #define SDL_IMPL
#define x264_IMPL

extern struct pstate pstate;
extern struct grid grid;
uint32_t *pixels;
bool quit;

void fill_pixels(uint32_t *pixels);

// test
#if 0
int main()
{
    int n = 32;
    init_sim(n);
    run(2.0);

    close_sim();
    return 0;
}

// principal impl
#else 


int main(int argc, char ** argv)
{   
    quit = false;

    #ifdef SDL_IMPL
        init_sdl();
    #else
        init_x264("out.mp4");
    #endif

    pixels = malloc(SCREEN_WIDTH * SCREEN_HEIGHT * sizeof(uint32_t));
    memset(pixels, 0xFF, SCREEN_HEIGHT * SCREEN_WIDTH * sizeof(uint32_t));
    ///

    u32 n = SCREEN_WIDTH;
    init_sim(n);
    //

    while (!quit)
    {
        if(grid.t >= 0.2)
            #ifdef SDL_IMPL
                init_problem(&pstate);
            #else
                break;
            #endif

        for (size_t i = 0; i < 10; i++)
            step(1.0);
        fill_pixels(pixels);

        //  for(int i = 0; i < SCREEN_WIDTH * SCREEN_HEIGHT; i++)
        //     pixels[i] = rand();

        #ifdef SDL_IMPL
            update_sdl();
        #else
            if(update_x264(pixels))
                break;
        #endif
    }

    #ifdef SDL_IMPL
        close_sdl();
    #else
        close_x264();
    #endif

    close_sim();
    free(pixels);

    return 0;
}

#endif

void fill_pixels(uint32_t *pixels)
{
    real *q = pstate.r;

    real min = q[grid.Ng], max = q[grid.Ng];
    // for (size_t i = grid.Ng; i < grid.Nx + grid.Ng; i++)
    // {
    //     // q[i] = (real)i / (real)grid.Nx;
    //     min = (min > q[i]) ? q[i] : min;
    //     max = (max < q[i]) ? q[i] : max;
    // }
    min = 0.0; max = 1.1;

    memset(pixels, 0xFF, SCREEN_WIDTH * SCREEN_HEIGHT * sizeof(uint32_t));
    size_t y = SCREEN_HEIGHT - 1;

    for (size_t i = 0; i < SCREEN_WIDTH; i++)
    {
        size_t old_y = y;
        y = (size_t)((1-(q[i]-min)/(max-min)) * SCREEN_HEIGHT);
        y = (y > SCREEN_HEIGHT - 1) ? SCREEN_HEIGHT - 1 : y;

        size_t min = (old_y < y) ? old_y : y;
        size_t max = (old_y > y) ? old_y : y;
        for(size_t j = min; j <= max; j++)
            pixels[i + j * SCREEN_WIDTH] = 0;
    }
}
