#include <stdio.h>
#include <stdlib.h>

#include "sim/sim.h"
#include "const.h"

#define SDL_IMPL
// #define x264_IMPL

#ifdef x264_IMPL
#include <signal.h>
#endif

extern struct pstate pstate;
extern struct grid grid;
uint32_t *pixels;
bool quit;

void fill_pixels(uint32_t *pixels);

// test

extern struct cstate cstate;
extern struct pstate x_slope;
extern struct pstate y_slope;

#if 0
int main()
{
    int n = 128;
    init_sim(n, n);
    // for(int i = 0;  i < 1000; i++){
    //     step(1.0);printf("%.18lf\n", grid.dt);
    // }
    // step(1.0);printf("%.18lf\n", grid.dt);
    // step(1.0);printf("%.18lf\n", grid.dt);
    // step(1.0);
    // printf("%.18lf\n", grid.dt);
    run(1.0);
    // conservative_to_primitive(&pstate, &cstate);

    // print_mat(pstate.r, n+4, n+4);
    // print_mat(x_slope.r, n+4, n+4);
    // step(1.0);
    // run(2.0);

    close_sim();
    return 0;
}

// principal impl
#else 

#ifdef x264_IMPL
void intHandler(int dummy) {
    close_x264();
    exit(0);
}

#endif

int main(int argc, char ** argv)
{   
    quit = false;

    #ifdef SDL_IMPL
        init_sdl();
    #else
        init_x264("out.mp4");
        signal(SIGINT, intHandler);
    #endif

    pixels = malloc(SCREEN_WIDTH * SCREEN_HEIGHT * sizeof(uint32_t));
    // memset(pixels, 0xFF, SCREEN_HEIGHT * SCREEN_WIDTH * sizeof(uint32_t));
    ///

    init_sim(SCREEN_WIDTH, SCREEN_WIDTH);
            // for (size_t i = 0; i < 20; i++)
            // {
            //     /* code */
            // step(1.0);
            // }
            
    //
    u32 itt = 0;
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
        // printf("%u\n", (++itt)*5);
        real dt = 1/200.;
        real tmax = grid.t + dt;
        while(grid.t < tmax)
            step(1.0);

        printf("%lf\n", grid.t);

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

u32 cmap_nipy_spectral(double v)
{
    static float data[22][3] = 
            {
                {  0.,   0.,   0.},
                {119.,   0., 136.},
                {136.,   0., 153.},
                {  0.,   0., 170.},
                {  0.,   0., 221.},
                {  0., 119., 221.},
                {  0., 153., 221.},
                {  0., 170., 170.},
                {  0., 170., 136.},
                {  0., 153.,   0.},
                {  0., 187.,   0.},
                {  0., 221.,   0.},
                {  0., 255.,   0.},
                {187., 255.,   0.},
                {238., 238.,   0.},
                {255., 204.,   0.},
                {255., 153.,   0.},
                {255.,   0.,   0.},
                {221.,   0.,   0.},
                {204.,   0.,   0.},
                {204., 204., 204.}
            };

    double index = v * 20.;
    double frac = index - floor(index);
    size_t id = index;

    u32 r = ((1-frac) * data[id][0] + frac * data[id+1][0]);
    u32 g = ((1-frac) * data[id][1] + frac * data[id+1][1]);
    u32 b = ((1-frac) * data[id][2] + frac * data[id+1][2]);

    return (r<<16) + (g<<8) + b;
}

void fill_pixels(uint32_t *pixels)
{
    real *q = pstate.r;
    double min = 0.8, max = 2.2;

    for(size_t j = 0; j < SCREEN_HEIGHT; j++)
        for (size_t i = 0; i < SCREEN_WIDTH; i++)
        {
            size_t id = (j+grid.Ng) * grid.Nx_tot + i + grid.Ng;
            
                // SOBEL
            #if 0
            int l = grid.Nx_tot;
            double vx = -   q[id]     +   q[id+2]
                        - 2*q[id+l]   + 2*q[id+l+2]
                        -   q[id+2*l] +   q[id+2*l+2];
            double vy = - q[id]     - 2*q[id+1]     - q[id+2]
                        + q[id+2*l] - 2*q[id+1+2*l] - q[id+2+2*l];
            double vv = sqrt(vx*vx+vy*vy);
            double v = (vv-min)/(max-min);
            uint8_t c = (v*256.);
            pixels[i + j * SCREEN_WIDTH] = (c<<16) + (c<<8) + c;

            #else
            double v = (q[id]-min)/(max-min);
            v = (v < 0.0) ? 0.0 : (v > 1.0) ? 1.0 : v;

            pixels[i + j * SCREEN_WIDTH] = cmap_nipy_spectral(1.-v);
            #endif
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
