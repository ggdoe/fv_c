#include "sim.h"

extern struct grid grid;
extern struct cstate cstate;

__attribute_maybe_unused__ 
void  fill_boundaries_step_periodic(real *q)
{
    // periodic
    u32 lo = grid.Ng;
    u32 hi = grid.Ng + grid.Nx;

    for(u32 i=0; i < grid.Ng; i++){
        q[i]      = q[hi - grid.Ng + i];
        q[hi + i] = q[lo + i];
    }
}

__attribute_maybe_unused__ 
void fill_boundaries_step_absorbing(real *q)
{
    // periodic
    u32 lo = grid.Ng;
    u32 hi = grid.Ng + grid.Nx;

    for(u32 i=0; i < grid.Ng; i++){
        q[i]      = q[lo];
        q[hi + i] = q[hi - grid.Ng + 1];
    }
}

__attribute_maybe_unused__ 
void fill_boundaries_step_reflect(real *q, double sign)
{
    // periodic.
    u32 lo = grid.Ng;
    u32 hi = grid.Ng + grid.Nx;

    for(u32 i=0; i < grid.Ng; i++){
        q[i]               = sign * q[lo + lo - 1 - i];
        q[hi - i + lo - 1] = sign * q[hi - lo + i];
    }
}

void fill_boundaries()
{
    // fill_boundaries_step_periodic(cstate.r);
    // fill_boundaries_step_periodic(cstate.ru);
    // fill_boundaries_step_periodic(cstate.e);

    fill_boundaries_step_absorbing(cstate.r);
    fill_boundaries_step_absorbing(cstate.ru);
    fill_boundaries_step_absorbing(cstate.e);

    // fill_boundaries_step_reflect(cstate.r ,  1.0);
    // fill_boundaries_step_reflect(cstate.ru, -1.0);
    // fill_boundaries_step_reflect(cstate.e ,  1.0);
}
