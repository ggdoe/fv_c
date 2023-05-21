#include <stdio.h>
#include <stdlib.h>

#include "sim/sim.h"
#include "const.h"

#define SDL_IMPL
// #define x264_IMPL

extern struct pstate pstate;
extern struct grid grid;
uint32_t *pixels;
bool quit;

void fill_pixels(uint32_t *pixels);

// test

extern struct cstate cstate;
extern struct pstate x_slope;
extern struct pstate y_slope;

#if 1
int main()
{
    int n = 8;
    init_sim(n, n);
    // for(int i = 0;  i < 200; i++)
    step(1.0);conservative_to_primitive(&pstate, &cstate);
    // printf("%.10lf\n", grid.dt);
    // step(1.0);conservative_to_primitive(&pstate, &cstate);
    // printf("%.10lf\n", grid.dt);
    // run(0.2);
    // conservative_to_primitive(&pstate, &cstate);
    print_mat(pstate.r, n+4, n+4);
    // print_mat(x_slope.r, n+4, n+4);
    // step(1.0);
    // run(2.0);

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

    init_sim(SCREEN_WIDTH, SCREEN_WIDTH);
            // for (size_t i = 0; i < 20; i++)
            // {
            //     /* code */
            // step(1.0);
            // }
            
    //

    while (!quit)
    {
            #if 1
                #ifdef SDL_IMPL
                    // if(grid.t >= 0.2)
                    //     init_problem(&pstate);
                    //     SDL_Delay(1000);
                #else
                    if(grid.t >= 1.0)
                        break;
                #endif
            #endif

        real dt = 1/60.;
        // real tmax = grid.t + dt;
        // while(grid.t < tmax)
        //     step(tmax - grid.t);
            step(1.0);
        fill_pixels(pixels);

        #ifdef SDL_IMPL
            update_sdl(pixels);
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
    size_t min = 0.0, max = 3.0;

    for(size_t j = 0; j < SCREEN_HEIGHT; j++)
        for (size_t i = 0; i < SCREEN_WIDTH; i++)
        {
            size_t id = (j) * grid.Nx_tot + i;
            int color = (int)(q[id]/max * 255);
            pixels[i + j * SCREEN_WIDTH] = (color<<16) + (color<<8) + color;
        }
}




// -----> old 1d //////

// void fill_pixels(uint32_t *pixels)
// {
//     real *q = pstate.r + SCREEN_WIDTH * SCREEN_HEIGHT / 2;

//     real min = q[grid.Ng], max = q[grid.Ng];
//     // for (size_t i = grid.Ng; i < grid.Nx + grid.Ng; i++)
//     // {
//     //     // q[i] = (real)i / (real)grid.Nx;
//     //     min = (min > q[i]) ? q[i] : min;
//     //     max = (max < q[i]) ? q[i] : max;
//     // }
//     min = 0.0; max = 1.1;

//     memset(pixels, 0xFF, SCREEN_WIDTH * SCREEN_HEIGHT * sizeof(uint32_t));
//     size_t y = (size_t)((1-(q[0]-min)/(max-min)) * SCREEN_HEIGHT);

//     for (size_t i = 0; i < SCREEN_WIDTH; i++)
//     {
//         size_t old_y = y;
//         y = (size_t)((1-(q[i]-min)/(max-min)) * SCREEN_HEIGHT);
//         y = (y > SCREEN_HEIGHT - 1) ? SCREEN_HEIGHT - 1 : y;

//         size_t min = (old_y < y) ? old_y : y;
//         size_t max = (old_y > y) ? old_y : y;
//         for(size_t j = min; j <= max; j++)
//             pixels[i + j * SCREEN_WIDTH] = 0;
//     }
// }
