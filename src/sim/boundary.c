#include "sim.h"

#ifdef _WIN32
#define __attribute_maybe_unused__
#endif

extern struct grid grid;

__attribute_maybe_unused__ 
void fill_boundaries_step_absorbing(real *q)
{
    // absorbing
    u32 Nx = grid.Nx_tot;
    u32 Ny = grid.Ny_tot;
    u32 lo = grid.Ng;
    u32 hix = Nx - grid.Ng;
    u32 hiy = Ny - grid.Ng;

    for (u32 y = 0; y < Ny; y++)
        for(u32 i=0; i < grid.Ng; i++)
        {
            q[Nx * y + i]       = q[Nx * y + lo];
            q[Nx * y + hix + i] = q[Nx * y + hix - 1];
        }

    for (u32 x = 0; x < Nx; x++)
        for(u32 j=0; j < grid.Ng; j++)
        {
            q[Nx * j + x]         = q[Nx * lo + x];
            q[Nx * (hiy + j) + x] = q[Nx * (hiy - 1) + x];
        }
}

__attribute_maybe_unused__ 
void  fill_boundaries_step_periodic(real *q)
{
    // periodic
    u32 Nx = grid.Nx_tot;
    u32 Ny = grid.Ny_tot;
    u32 lo = grid.Ng;
    u32 hix = Nx - grid.Ng;
    u32 hiy = Ny - grid.Ng;

    for (u32 y = 0; y < Ny; y++)
        for(u32 i=0; i < grid.Ng; i++)
        {
            q[Nx * y + i]       = q[Nx * y + hix - lo + i];
            q[Nx * y + hix + i] = q[Nx * y + lo + i];
        }

    for (u32 x = 0; x < Nx; x++)
        for(u32 j=0; j < grid.Ng; j++)
        {
            q[Nx * j + x]         = q[Nx * (hiy - lo + j) + x];
            q[Nx * (hiy + j) + x] = q[Nx * (lo + j) + x];
        }
}

__attribute_maybe_unused__ 
void fill_boundaries_step_reflect(real *q, double sign)
{
    // reflect
    u32 Nx = grid.Nx_tot;
    u32 Ny = grid.Ny_tot;
    u32 lo = grid.Ng;
    u32 hix = Nx - grid.Ng;
    u32 hiy = Ny - grid.Ng;

    for (u32 y = 0; y < Ny; y++)
        for(u32 i=0; i < grid.Ng; i++)
        {
            q[Nx * y + i]       = sign * q[Nx * y + lo + lo - i - 1];
            q[Nx * y + hix + i] = sign * q[Nx * y + hix - 1 - i];
        }

    for (u32 x = 0; x < Nx; x++)
        for(u32 j=0; j < grid.Ng; j++)
        {
            q[Nx * j + x]         = sign * q[Nx * (lo + lo - j - 1) + x];
            q[Nx * (hiy + j) + x] = sign * q[Nx * (hiy - 1 - j) + x];
        }
}

void fill_boundaries(struct cstate *cstate)
{
    // fill_boundaries_step_periodic(cstate->r);
    // fill_boundaries_step_periodic(cstate->ru);
    // fill_boundaries_step_periodic(cstate->rv);
    // fill_boundaries_step_periodic(cstate->e);

    // fill_boundaries_step_absorbing(cstate->r);
    // fill_boundaries_step_absorbing(cstate->ru);
    // fill_boundaries_step_absorbing(cstate->rv);
    // fill_boundaries_step_absorbing(cstate->e);

    fill_boundaries_step_reflect(cstate->r ,  1.0);
    fill_boundaries_step_reflect(cstate->ru, -1.0);
    fill_boundaries_step_reflect(cstate->rv, -1.0);
    fill_boundaries_step_reflect(cstate->e ,  1.0);
}
