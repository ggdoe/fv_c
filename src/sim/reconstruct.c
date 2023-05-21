#include "sim.h"

extern struct grid grid;
extern struct pstate pstate;
struct pstate x_slope;

real minmod(real a, real b)
{
    return (a * b < 0.0) ? 0.0 : (fabs(a) < fabs(b)) ? a : b;
}

void compute_slopes()
{
    for (size_t i = 1; i < grid.Nx_tot - 1; i++)
    {
        x_slope.r[i] = minmod(pstate.r[i  ] - pstate.r[i-1], 
                              pstate.r[i+1] - pstate.r[i  ]);
        x_slope.u[i] = minmod(pstate.u[i  ] - pstate.u[i-1], 
                              pstate.u[i+1] - pstate.u[i  ]);
        x_slope.p[i] = minmod(pstate.p[i  ] - pstate.p[i-1], 
                              pstate.p[i+1] - pstate.p[i  ]);
    }

    x_slope.r[0]               = x_slope.r[1]; 
    x_slope.u[0]               = x_slope.u[1]; 
    x_slope.p[0]               = x_slope.p[1]; 
    x_slope.r[grid.Nx_tot - 1] = x_slope.r[grid.Nx_tot - 2];
    x_slope.u[grid.Nx_tot - 1] = x_slope.u[grid.Nx_tot - 2];
    x_slope.p[grid.Nx_tot - 1] = x_slope.p[grid.Nx_tot - 2];
}

void reconstruct_interface(struct pcell *p, size_t i, real sign)
{
        p->r = p->r + sign * 0.5 * x_slope.r[i];
        p->u = p->u + sign * 0.5 * x_slope.u[i];
        p->p = p->p + sign * 0.5 * x_slope.p[i];
}


void muscl_hancock(struct pcell *p, size_t i, real sign)
{
    // TODO après 2d
    
    // p->r = p->r + sign * 0.5 * x_slope.r[i];
    // p->u = p->u + sign * 0.5 * x_slope.u[i];
    // p->p = p->p + sign * 0.5 * x_slope.p[i];
}
