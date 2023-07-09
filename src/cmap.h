#include "math.h"
#include "const.h"

#include "sim/sim.h"

extern struct grid grid;

static inline
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
                {204., 204., 204.},
                {204., 204., 204.}
            };

    v = (v < 0.0) ? 0.0 : (v > 1.0) ? 1.0 : v;
    double index = v * 20.;
    double frac = index - floor(index);
    size_t id = index;

    u32 r = ((1-frac) * data[id][0] + frac * data[id+1][0]);
    u32 g = ((1-frac) * data[id][1] + frac * data[id+1][1]);
    u32 b = ((1-frac) * data[id][2] + frac * data[id+1][2]);

    return (r<<16) | (g<<8) | b;
}

static inline
u32 sobel(real *q, size_t id, double min, double max)
{
        int l = grid.Nx_tot;
        double vx = -   q[id]     +   q[id+2]
                    - 2*q[id+l]   + 2*q[id+l+2]
                    -   q[id+2*l] +   q[id+2*l+2];
        double vy = - q[id]     - 2*q[id+1]     - q[id+2]
                    + q[id+2*l] - 2*q[id+1+2*l] - q[id+2+2*l];
        double d = sqrt(vx*vx+vy*vy);
        double v = (d-min)/(max-min);
        uint8_t c = (uint8_t)(v*255.);
        return (c<<16) | (c<<8) | c;
}

void fill_pixels(uint32_t *pixels, real *q)
{
    // double min = 0.8, max = 2.2;
    double min = 0.08, max = 6.5;

    for(size_t j = 0; j < SCREEN_HEIGHT; j++)
        for (size_t i = 0; i < SCREEN_WIDTH; i++)
        {
            size_t id = (j+grid.Ng) * grid.Nx_tot + i + grid.Ng;
            
                // SOBEL
            #if 0
            pixels[i + j * SCREEN_WIDTH] = sobel(q, id, min, max);
            #else
            double v = (q[id]-min)/(max-min);
            pixels[i + j * SCREEN_WIDTH] = cmap_nipy_spectral(v);
            #endif
        }
}
