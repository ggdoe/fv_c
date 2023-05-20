#include "sim.h"

struct grid grid;

struct pstate pstate;
struct cstate cstate;
struct fluxes fluxes;

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

    // x_slope.r[0]               = x_slope.r[1]; 
    // x_slope.u[0]               = x_slope.u[1]; 
    // x_slope.p[0]               = x_slope.p[1]; 
    // x_slope.r[grid.Nx_tot - 1] = x_slope.r[grid.Nx_tot - 2];
    // x_slope.u[grid.Nx_tot - 1] = x_slope.u[grid.Nx_tot - 2];
    // x_slope.p[grid.Nx_tot - 1] = x_slope.p[grid.Nx_tot - 2];
}

void reconstruct_interface(struct pcell *p, size_t i, real sign)
{
        p->r = p->r + sign * 0.5 * x_slope.r[i];
        p->u = p->u + sign * 0.5 * x_slope.u[i];
        p->p = p->p + sign * 0.5 * x_slope.p[i];
}

void update_cells()
{
    for (size_t i = 0; i < grid.Nx_tot - 1; i++)
    {
        struct pcell pl = get_cell(&pstate, i);
        struct pcell pr = get_cell(&pstate, i+1);
        reconstruct_interface(&pl, i,   +1.0);
        reconstruct_interface(&pr, i+1, -1.0);

        // printf("%lu\n",i);
        struct fcell F  = solve_fluxes(&pl, &pr);
        
        fluxes.r[i]  = F.r;
        fluxes.ru[i] = F.ru;
        fluxes.e[i]  = F.e;
    }

    for (size_t i = 0; i < grid.Nx_tot-2; i++)
    {
        real dtdx = grid.dt / grid.dx;
        cstate.r[i+1]  += dtdx * (fluxes.r[i]  - fluxes.r[i+1]);
        cstate.ru[i+1] += dtdx * (fluxes.ru[i] - fluxes.ru[i+1]);
        cstate.e[i+1]  += dtdx * (fluxes.e[i]  - fluxes.e[i+1]);
    }
}

void run(real tmax)
{
    real *t = &grid.t;
    // size_t itt = 0;

    while(*t < tmax)
    {
        step(tmax-*t);
        // itt++;
        // printf("%lu - %lf\n", itt, t);
    }

}

void compute_dt()
{
    real dt = __FLT_MAX__;
    for (size_t i = 0; i < grid.Nx_tot; i++)
    {
        real cs = SQRT(pstate.p[i] * grid.gamma / pstate.r[i]);
        real new_dt = grid.CFL * grid.dx / (fabs(cs) + fabs(pstate.u[i]));
        dt = (new_dt < dt) ? new_dt : dt;
    }
    // printf("%lf\n", dt);
    grid.dt = dt;
}

void step(real dt_max)
{
    compute_dt();
    if(grid.dt > dt_max)
        grid.dt = dt_max;
    grid.t += grid.dt;

    compute_slopes();
    update_cells();
    fill_boundaries();
    conservative_to_primitive(&pstate, &cstate);
}
