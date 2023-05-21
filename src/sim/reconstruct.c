#include "sim.h"

extern struct grid grid;
extern struct pstate pstate;
struct pstate x_slope;
struct pstate y_slope;

real minmod(real a, real b)
{
    return (a * b < 0.0) ? 0.0 : (fabs(a) < fabs(b)) ? a : b;
}

void compute_slopes()
{
    for (size_t j = 1; j < grid.Ny_tot - 1; j++)
        for (size_t i = 1; i < grid.Nx_tot - 1; i++)
        {
            size_t id = j * grid.Nx_tot + i;
            x_slope.r[id] = minmod(pstate.r[id  ] - pstate.r[id-1], 
                                   pstate.r[id+1] - pstate.r[id  ]);
            x_slope.u[id] = minmod(pstate.u[id  ] - pstate.u[id-1], 
                                   pstate.u[id+1] - pstate.u[id  ]);
            x_slope.p[id] = minmod(pstate.p[id  ] - pstate.p[id-1], 
                                   pstate.p[id+1] - pstate.p[id  ]);
                                   
            size_t l = grid.Nx_tot; // access top/bot cell
            y_slope.r[id] = minmod(pstate.r[id  ] - pstate.r[id-l], 
                                   pstate.r[id+l] - pstate.r[id  ]);
            y_slope.u[id] = minmod(pstate.u[id  ] - pstate.u[id-l], 
                                   pstate.u[id+l] - pstate.u[id  ]);
            y_slope.p[id] = minmod(pstate.p[id  ] - pstate.p[id-l], 
                                   pstate.p[id+l] - pstate.p[id  ]);
        }

    // x_slope.r[0]               = x_slope.r[1]; 
    // x_slope.u[0]               = x_slope.u[1]; 
    // x_slope.p[0]               = x_slope.p[1]; 
    // x_slope.r[grid.Nx_tot - 1] = x_slope.r[grid.Nx_tot - 2];
    // x_slope.u[grid.Nx_tot - 1] = x_slope.u[grid.Nx_tot - 2];
    // x_slope.p[grid.Nx_tot - 1] = x_slope.p[grid.Nx_tot - 2];
}

void reconstruct_interface_x(struct pcell *p, size_t i, real sign)
{
        p->r = p->r + sign * 0.5 * x_slope.r[i];
        p->u = p->u + sign * 0.5 * x_slope.u[i];
        p->p = p->p + sign * 0.5 * x_slope.p[i];
}

void reconstruct_interface_y(struct pcell *p, size_t i, real sign)
{
        p->r = p->r + sign * 0.5 * y_slope.r[i];
        p->u = p->u + sign * 0.5 * y_slope.u[i];
        p->p = p->p + sign * 0.5 * y_slope.p[i];
}


void muscl_hancock(struct pcell *p, size_t i, real sign)
{
    // TODO après 2d
    
    // p->r = p->r + sign * 0.5 * x_slope.r[i];
    // p->u = p->u + sign * 0.5 * x_slope.u[i];
    // p->p = p->p + sign * 0.5 * x_slope.p[i];
}
